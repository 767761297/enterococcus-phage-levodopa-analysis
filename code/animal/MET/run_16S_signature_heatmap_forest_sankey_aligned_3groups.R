#!/usr/bin/env Rscript

# ============================================================
# Sankey-aligned 16S fingerprint workflow (3 groups)
# Combined script for:
#   1) Signature / fingerprint taxa identification (6 taxonomic levels)
#   2) Fingerprint heatmaps (6 taxonomic levels)
#   3) Forest plots for fingerprint taxa (6 taxonomic levels; GR vs MR)
#
# Statistical alignment target:
#   species–metabolite Sankey species-side workflow
#   - prevalence filter across retained 3-group samples
#   - mean relative abundance filter across retained 3-group samples
#   - CLR on counts with pseudocount = 0.5
#   - Wilcoxon rank-sum test
#   - BH-FDR
#
# Comparison scope remains unchanged from legacy scripts:
#   - signature / forest testing: GR vs MR only
#   - heatmap display: GR / MR / MR+Phage cocktail
#
# Input:
#   allsample_bundle_from_featuretable_3groups/allsample.<level>.abundance.xls
#
# Organized outputs:
#   16S_SignatureHeatmapForest_sankey_aligned_3groups/
#     01_signature_q0.1_6levels/
#     02_heatmaps_q0.1_6levels/
#     03_forest_q0.1_6levels/
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(grid)
})

# -------------------- USER SETTINGS --------------------
IN_DIR <- "allsample_bundle_from_featuretable_3groups"
ROOT_OUTDIR <- "16S_SignatureHeatmapForest_sankey_aligned_3groups"

SIG_OUTDIR    <- file.path(ROOT_OUTDIR, "01_signature_q0.1_6levels")
HEAT_OUTDIR   <- file.path(ROOT_OUTDIR, "02_heatmaps_q0.1_6levels")
FOREST_OUTDIR <- file.path(ROOT_OUTDIR, "03_forest_q0.1_6levels")

for (d in c(ROOT_OUTDIR, SIG_OUTDIR, HEAT_OUTDIR, FOREST_OUTDIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

LEVELS <- c("phylum", "class", "order", "family", "genus", "species")
LEVELS_DISPLAY <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
names(LEVELS_DISPLAY) <- LEVELS

# three-group sample order
SAMPLE_ORDER <- c(paste0("A", 1:5), paste0("B", 1:5), paste0("C", 1:5))
GROUP_LEVELS <- c("GR", "MR", "MR+Phage cocktail")

map_group <- function(s) {
  if (str_detect(s, "^A")) return("GR")
  if (str_detect(s, "^B")) return("MR")
  if (str_detect(s, "^C")) return("MR+Phage cocktail")
  return(NA_character_)
}

GROUP_COLORS <- c(
  "GR" = "#7CAE00",
  "MR" = "#F8766D",
  "MR+Phage cocktail" = "#00BFC4"
)

# Sankey-aligned statistics
PSEUDO <- 0.5
PREV_CUT <- 0.10
MEAN_RELAB_CUT <- 1e-4
Q_CUT <- 0.10

# Heatmap options
CLUSTER_COLS <- FALSE
SHOW_COLNAMES <- FALSE
FONTSIZE_ROW <- 7

# Forest options
USE_SIGNATURE_FILES <- TRUE
BOOT_B <- 800
SEED <- 20260106
COL_UP   <- "#D7191C"   # GR > MR
COL_DOWN <- "#2C7BB6"   # GR < MR
# ------------------------------------------------------

# -------------------- HELPERS --------------------
log_msg <- function(...) cat(paste0(..., "\n"))

read_abundance_table <- function(level) {
  f1 <- file.path(IN_DIR, paste0("allsample.", level, ".abundance.xls"))
  f2 <- paste0("allsample.", level, ".abundance.xls")
  f <- if (file.exists(f1)) f1 else f2
  if (!file.exists(f)) stop("Cannot find input file for level '", level, "': ", f1, " or ", f2)

  df <- read.delim(f, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(ncol(df) >= 2)

  tax_col <- colnames(df)[1]
  df <- df %>% rename(taxon = all_of(tax_col))

  sample_cols <- intersect(colnames(df), SAMPLE_ORDER)
  if (length(sample_cols) == 0) {
    sample_cols <- setdiff(colnames(df), "taxon")
    sample_cols <- sample_cols[grepl("^[A-C][0-9]+$", sample_cols)]
  }
  if (length(sample_cols) == 0) stop("No sample columns detected in: ", f)

  mat <- as.matrix(df[, sample_cols, drop = FALSE])
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  rownames(mat) <- make.unique(as.character(df$taxon))

  list(mat = mat, file = f)
}

order_and_keep_3groups <- function(counts) {
  samples <- colnames(counts)
  meta <- tibble(sample = samples, group = vapply(samples, map_group, character(1))) %>%
    filter(group %in% GROUP_LEVELS)

  keep_samp <- meta$sample
  counts <- counts[, keep_samp, drop = FALSE]

  samp2 <- SAMPLE_ORDER[SAMPLE_ORDER %in% colnames(counts)]
  counts <- counts[, samp2, drop = FALSE]
  meta <- meta %>%
    mutate(group = factor(group, levels = GROUP_LEVELS)) %>%
    arrange(match(sample, colnames(counts)))

  list(counts = counts, meta = meta)
}

filter_counts_sankey_aligned <- function(counts) {
  if (nrow(counts) == 0) return(counts)
  rel <- sweep(counts, 2, colSums(counts), "/")
  rel[!is.finite(rel)] <- 0

  prev <- rowMeans(counts > 0)
  mean_relab <- rowMeans(rel)
  keep <- prev >= PREV_CUT & mean_relab >= MEAN_RELAB_CUT
  counts[keep, , drop = FALSE]
}

clr_transform <- function(counts, pseudo = PSEUDO) {
  x <- counts + pseudo
  logx <- log(x)
  sweep(logx, 2, colMeans(logx), "-")
}

safe_wilcox <- function(x1, x0) {
  idx <- is.finite(x1) & is.finite(x0)
  if (length(x1) < 2 || length(x0) < 2) return(NA_real_)
  if (length(unique(c(x1, x0))) <= 1) return(1)
  out <- tryCatch(
    suppressWarnings(wilcox.test(x1, x0, exact = FALSE)$p.value),
    error = function(e) NA_real_
  )
  out
}

calc_log2fc_rel <- function(counts, s1, s0, eps = 1e-6) {
  rel <- sweep(counts, 2, colSums(counts), "/")
  rel[!is.finite(rel)] <- 0
  m1 <- rowMeans(rel[, s1, drop = FALSE]) + eps
  m0 <- rowMeans(rel[, s0, drop = FALSE]) + eps
  log2(m1 / m0)
}

bootstrap_ci <- function(x_gr, x_mr, B = 800, seed = 1) {
  set.seed(seed)
  ng <- length(x_gr)
  nm <- length(x_mr)
  if (ng < 2 || nm < 2) return(c(NA_real_, NA_real_))
  boot <- replicate(B, {
    mean(sample(x_gr, ng, replace = TRUE)) - mean(sample(x_mr, nm, replace = TRUE))
  })
  stats::quantile(boot, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
}

read_signature_taxa <- function(level) {
  f <- file.path(SIG_OUTDIR, paste0("SIGNATURE_q0.1_", level, "_GR_vs_MR.csv"))
  if (!file.exists(f)) return(character(0))
  sig <- readr::read_csv(f, show_col_types = FALSE)
  if (!("taxon" %in% names(sig))) return(character(0))
  unique(as.character(sig$taxon))
}

# -------------------- STEP 1: SIGNATURE --------------------
log_msg("[STEP1] Build 16S signatures (Sankey-aligned stats)")
summary_rows <- list()

for (lv in LEVELS) {
  log_msg("\n[RUN] level = ", lv)

  obj <- read_abundance_table(lv)
  tmp <- order_and_keep_3groups(obj$mat)
  counts <- tmp$counts
  meta <- tmp$meta

  counts <- filter_counts_sankey_aligned(counts)
  if (nrow(counts) == 0) {
    log_msg("[SKIP] ", lv, ": no taxa remain after Sankey-aligned filtering.")
    next
  }

  s_GR <- meta$sample[meta$group == "GR"]
  s_MR <- meta$sample[meta$group == "MR"]
  if (length(s_GR) < 2 || length(s_MR) < 2) stop("Not enough samples in GR or MR for: ", lv)

  clr <- clr_transform(counts, pseudo = PSEUDO)

  delta_clr <- rowMeans(clr[, s_GR, drop = FALSE]) - rowMeans(clr[, s_MR, drop = FALSE])
  pvals <- apply(clr, 1, function(v) safe_wilcox(v[s_GR], v[s_MR]))
  qvals <- p.adjust(pvals, method = "BH")
  log2fc_rel <- calc_log2fc_rel(counts, s_GR, s_MR)

  rel <- sweep(counts, 2, colSums(counts), "/")
  rel[!is.finite(rel)] <- 0
  mean_relab_GR <- rowMeans(rel[, s_GR, drop = FALSE])
  mean_relab_MR <- rowMeans(rel[, s_MR, drop = FALSE])
  mean_relab_all <- rowMeans(rel)

  prev_all <- rowMeans(counts > 0)
  prev_GR  <- rowMeans(counts[, s_GR, drop = FALSE] > 0)
  prev_MR  <- rowMeans(counts[, s_MR, drop = FALSE] > 0)

  all_tbl <- tibble(
    level = lv,
    taxon = rownames(counts),
    delta_CLR = as.numeric(delta_clr),
    log2FC_rel = as.numeric(log2fc_rel),
    mean_relab_all = as.numeric(mean_relab_all),
    mean_relab_GR = as.numeric(mean_relab_GR),
    mean_relab_MR = as.numeric(mean_relab_MR),
    prev_all = as.numeric(prev_all),
    prev_GR  = as.numeric(prev_GR),
    prev_MR  = as.numeric(prev_MR),
    p_value = as.numeric(pvals),
    q_value = as.numeric(qvals),
    direction = case_when(
      delta_CLR > 0 ~ "Higher_in_GR",
      delta_CLR < 0 ~ "Lower_in_GR",
      TRUE ~ "Zero"
    )
  ) %>%
    arrange(q_value, desc(abs(delta_CLR)))

  out_all <- file.path(SIG_OUTDIR, paste0("ALL_", lv, "_GR_vs_MR_CLR_Wilcox_BH.csv"))
  write_csv(all_tbl, out_all)

  sig_tbl <- all_tbl %>% filter(!is.na(q_value), q_value < Q_CUT)
  out_sig <- file.path(SIG_OUTDIR, paste0("SIGNATURE_q0.1_", lv, "_GR_vs_MR.csv"))
  write_csv(sig_tbl, out_sig)

  summary_rows[[lv]] <- tibble(
    level = lv,
    input_file = obj$file,
    n_tested = nrow(all_tbl),
    n_signature_q0.1 = nrow(sig_tbl),
    q_min = suppressWarnings(min(all_tbl$q_value, na.rm = TRUE))
  )

  log_msg(sprintf("[OK] %s: tested=%d, signature(q<0.1)=%d, q_min=%.4g",
                  lv, nrow(all_tbl), nrow(sig_tbl),
                  suppressWarnings(min(all_tbl$q_value, na.rm = TRUE))))
}

summary_tbl <- bind_rows(summary_rows)
write_csv(summary_tbl, file.path(SIG_OUTDIR, "SUMMARY_signature_q0.1_6levels.csv"))
log_msg("\n[DONE] STEP1 outputs in: ", SIG_OUTDIR)

# -------------------- STEP 2: HEATMAP --------------------
log_msg("\n[STEP2] Build 6-level fingerprint heatmaps")

make_heatmap_one_level <- function(level) {
  sig_taxa <- read_signature_taxa(level)
  if (length(sig_taxa) == 0) {
    message("[SKIP] ", level, ": signature set is empty (q<0.1 produced 0 taxa).")
    return(invisible(NULL))
  }

  obj <- read_abundance_table(level)
  tmp <- order_and_keep_3groups(obj$mat)
  counts <- tmp$counts
  meta <- tmp$meta

  counts <- filter_counts_sankey_aligned(counts)
  keep <- intersect(sig_taxa, rownames(counts))
  if (length(keep) == 0) {
    message("[SKIP] ", level, ": none of signature taxa found after Sankey-aligned filtering.")
    return(invisible(NULL))
  }

  counts_sig <- counts[keep, , drop = FALSE]
  clr_sig <- clr_transform(counts_sig, PSEUDO)

  samples <- colnames(clr_sig)
  groups  <- vapply(samples, map_group, character(1))
  groups  <- factor(groups, levels = GROUP_LEVELS)

  anno_col <- data.frame(group = groups)
  rownames(anno_col) <- samples

  present_groups <- names(table(groups))
  gaps_col <- cumsum(as.integer(table(groups)))
  if (length(gaps_col) > 0) gaps_col <- gaps_col[-length(gaps_col)]

  anno_colors <- list(group = GROUP_COLORS[present_groups])

  row_sd <- apply(clr_sig, 1, sd)
  clr_sig <- clr_sig[row_sd > 0, , drop = FALSE]
  if (nrow(clr_sig) == 0) {
    message("[SKIP] ", level, ": all signature taxa have zero variance after filtering.")
    return(invisible(NULL))
  }

  h_in <- max(4, min(14, 0.22 * nrow(clr_sig) + 2.5))
  w_in <- 8.8

  out_png <- file.path(HEAT_OUTDIR, paste0("Heatmap_signature_q0.1_", level, "_3groups.png"))
  out_pdf <- file.path(HEAT_OUTDIR, paste0("Heatmap_signature_q0.1_", level, "_3groups.pdf"))

  pheatmap(
    clr_sig,
    scale = "row",
    annotation_col = anno_col,
    annotation_colors = anno_colors,
    cluster_rows = TRUE,
    cluster_cols = CLUSTER_COLS,
    gaps_col = if (!CLUSTER_COLS && length(gaps_col) > 0) gaps_col else NULL,
    show_colnames = SHOW_COLNAMES,
    fontsize_row = FONTSIZE_ROW,
    fontsize_col = 8,
    main = NA,
    filename = out_png,
    width = w_in, height = h_in
  )

  pheatmap(
    clr_sig,
    scale = "row",
    annotation_col = anno_col,
    annotation_colors = anno_colors,
    cluster_rows = TRUE,
    cluster_cols = CLUSTER_COLS,
    gaps_col = if (!CLUSTER_COLS && length(gaps_col) > 0) gaps_col else NULL,
    show_colnames = SHOW_COLNAMES,
    fontsize_row = FONTSIZE_ROW,
    fontsize_col = 8,
    main = NA,
    filename = out_pdf,
    width = w_in, height = h_in
  )

  message(sprintf("[OK] %s: taxa=%d -> %s", level, nrow(clr_sig), out_png))
  invisible(list(level = level, n_taxa = nrow(clr_sig), png = out_png, pdf = out_pdf))
}

res_heat <- lapply(LEVELS, make_heatmap_one_level)

merged_pdf <- file.path(HEAT_OUTDIR, "Heatmap_signature_q0.1_6levels_merged_3groups.pdf")
pdf(merged_pdf, width = 9.5, height = 8)
for (lv in LEVELS) {
  sig_taxa <- tryCatch(read_signature_taxa(lv), error = function(e) character(0))
  if (length(sig_taxa) == 0) next

  obj <- tryCatch(read_abundance_table(lv), error = function(e) NULL)
  if (is.null(obj)) next

  tmp <- order_and_keep_3groups(obj$mat)
  counts <- filter_counts_sankey_aligned(tmp$counts)

  keep <- intersect(sig_taxa, rownames(counts))
  if (length(keep) == 0) next

  clr_sig <- clr_transform(counts[keep, , drop = FALSE], PSEUDO)
  row_sd <- apply(clr_sig, 1, sd)
  clr_sig <- clr_sig[row_sd > 0, , drop = FALSE]
  if (nrow(clr_sig) == 0) next

  samples <- colnames(clr_sig)
  groups  <- factor(vapply(samples, map_group, character(1)), levels = GROUP_LEVELS)
  anno_col <- data.frame(group = groups)
  rownames(anno_col) <- samples

  present_groups <- names(table(groups))
  gaps_col <- cumsum(as.integer(table(groups)))
  if (length(gaps_col) > 0) gaps_col <- gaps_col[-length(gaps_col)]

  anno_colors <- list(group = GROUP_COLORS[present_groups])

  ph <- pheatmap(
    clr_sig,
    scale = "row",
    annotation_col = anno_col,
    annotation_colors = anno_colors,
    cluster_rows = TRUE,
    cluster_cols = CLUSTER_COLS,
    gaps_col = if (!CLUSTER_COLS && length(gaps_col) > 0) gaps_col else NULL,
    show_colnames = SHOW_COLNAMES,
    fontsize_row = FONTSIZE_ROW,
    fontsize_col = 8,
    main = NA
  )
  grid.newpage()
  grid.draw(ph$gtable)
}
dev.off()

log_msg("[DONE] STEP2 heatmaps saved in: ", HEAT_OUTDIR)
log_msg("[DONE] STEP2 merged PDF: ", merged_pdf)

# -------------------- STEP 3: FOREST --------------------
log_msg("\n[STEP3] Forest plots for fingerprint taxa (GR vs MR)")
all_tables <- list()

for (lv in LEVELS) {
  lv_disp <- LEVELS_DISPLAY[[lv]]
  log_msg("\n[LEVEL] ", lv_disp, " ...")

  obj <- read_abundance_table(lv)
  tmp <- order_and_keep_3groups(obj$mat)
  counts <- filter_counts_sankey_aligned(tmp$counts)
  meta <- tmp$meta

  s_GR <- meta$sample[meta$group == "GR"]
  s_MR <- meta$sample[meta$group == "MR"]
  if (length(s_GR) < 2 || length(s_MR) < 2) stop("Need >=2 samples in both GR and MR.")

  if (nrow(counts) < 3) {
    log_msg("[SKIP] too few taxa after Sankey-aligned filtering at ", lv_disp)
    next
  }

  clr <- clr_transform(counts, pseudo = PSEUDO)

  taxa <- rownames(clr)
  pvals <- vapply(taxa, function(t) {
    xg <- clr[t, s_GR]
    xm <- clr[t, s_MR]
    safe_wilcox(xg, xm)
  }, numeric(1))
  qvals <- p.adjust(pvals, method = "BH")
  eff <- vapply(taxa, function(t) mean(clr[t, s_GR]) - mean(clr[t, s_MR]), numeric(1))

  stats_all <- data.frame(
    level = lv_disp,
    taxon = taxa,
    effect = eff,
    p_value = pvals,
    q_value = qvals,
    direction = ifelse(eff >= 0, "GR higher", "MR higher"),
    stringsAsFactors = FALSE
  )

  sig_taxa <- character(0)
  if (USE_SIGNATURE_FILES) sig_taxa <- read_signature_taxa(lv)

  if (length(sig_taxa) > 0) {
    sig <- stats_all %>% filter(taxon %in% sig_taxa)
    log_msg("[INFO] signature from file: n=", nrow(sig))
  } else {
    sig <- stats_all %>% filter(q_value < Q_CUT)
    log_msg("[INFO] signature from q<", Q_CUT, ": n=", nrow(sig))
  }

  if (nrow(sig) == 0) {
    log_msg("[SKIP] no signature taxa at ", lv_disp)
    next
  }

  set.seed(SEED)
  cis <- t(vapply(sig$taxon, function(t) {
    xg <- clr[t, s_GR]
    xm <- clr[t, s_MR]
    bootstrap_ci(xg, xm, B = BOOT_B, seed = SEED + match(t, sig$taxon))
  }, numeric(2)))
  sig$ci_low  <- cis[, 1]
  sig$ci_high <- cis[, 2]

  sig <- sig %>%
    arrange(effect) %>%
    mutate(taxon = factor(taxon, levels = taxon))

  out_csv <- file.path(FOREST_OUTDIR, paste0("Forest_signature_q0.1_", lv, "_GR_vs_MR.csv"))
  write.csv(sig, out_csv, row.names = FALSE)

  all_tables[[lv_disp]] <- sig

  n_tax <- nrow(sig)
  h <- max(4.5, 0.22 * n_tax + 1.2)

  p <- ggplot(sig, aes(x = effect, y = taxon)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey40") +
    geom_errorbar(aes(xmin = ci_low, xmax = ci_high, color = direction), height = 0.2, linewidth = 0.7, orientation = "y") +
    geom_point(aes(color = direction), size = 2.4, shape = 16) +
    scale_color_manual(values = c("GR higher" = COL_UP, "MR higher" = COL_DOWN)) +
    theme_classic(base_size = 12) +
    theme(
      legend.title = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y  = element_text(size = 9)
    ) +
    labs(x = "Effect size (mean CLR difference)")

  out_pdf <- file.path(FOREST_OUTDIR, paste0("Forest_signature_q0.1_", lv, "_GR_vs_MR.pdf"))
  out_png <- file.path(FOREST_OUTDIR, paste0("Forest_signature_q0.1_", lv, "_GR_vs_MR.png"))

  ggsave(out_pdf, p, width = 8.2, height = h)
  ggsave(out_png, p, width = 8.2, height = h, dpi = 300)

  log_msg("[OK] ", lv_disp, " -> ", basename(out_pdf), " | table: ", basename(out_csv))
}

if (length(all_tables) > 0 && requireNamespace("writexl", quietly = TRUE)) {
  out_xlsx <- file.path(FOREST_OUTDIR, "Forest_signature_q0.1_6Levels_GR_vs_MR.xlsx")
  writexl::write_xlsx(all_tables, out_xlsx)
  log_msg("\n[OK] Combined xlsx: ", out_xlsx)
} else {
  log_msg("\n[NOTE] writexl not installed -> skip combined xlsx (CSV per level already written).")
}

log_msg("\n[DONE] STEP3 outputs in: ", FOREST_OUTDIR)
log_msg("\n[ALL DONE] Root output dir: ", ROOT_OUTDIR)
