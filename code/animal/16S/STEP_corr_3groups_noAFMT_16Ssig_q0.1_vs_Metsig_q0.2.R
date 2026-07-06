## ============================================================
## STEP_corr_3groups_noAFMT_16Ssig_q0.1_vs_Metsig_q0.2.R
## Keep 3 groups only: GR (A*), MR (B*), MR+Phage cocktail (C*)
## Drop MR+AFMT (D*) from ALL steps (16S fingerprint + correlation).
##
## Met fingerprint: from neurotransmitter_metabolomics_stats_q0_2.xlsx (all_results),
##                  q<0.2 for GR vs MR.
## Met quant: Q2.xlsx, Compound Name (nmol/g) + A1..D5 columns
## 16S: search allsample.<level>.reabundance/abundance under . and allsample_bundle_from_featuretable
##
## Outputs:
## - MetFingerprint_GRvsMR_q0.20.csv
## - 16S_<level>_GRvsMR_Fingerprint_q0.10_3groups_noAFMT.csv
## - CORR_<level>_..._3groups_noAFMT.csv  (overall + group_adjusted)
## - EDGES_<level>_..._3groups_noAFMT.csv
## ============================================================

options(error = NULL)

## ---------- packages ----------
need_pkgs <- c("readxl", "dplyr", "stringr", "tidyr", "purrr", "readr")
for (p in need_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop("Missing package: ", p, ". Install with install.packages('", p, "')")
  }
}
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(readr)
})
HAS_PHEATMAP <- requireNamespace("pheatmap", quietly = TRUE)

## ---------- user settings ----------
MET_QUANT_FILE  <- "Q2.xlsx"
MET_STATS_FILE  <- "neurotransmitter_metabolomics_stats_q0_2.xlsx"
MET_STATS_SHEET <- "all_results"

MET_Q_CUT <- 0.2    # metabolite fingerprint threshold (GR vs MR)
TAX_Q_CUT <- 0.1    # 16S fingerprint threshold (GR vs MR)

METHOD_COR <- "spearman"

# correlation filtering for “edges/heatmap”
CORR_Q_CUT   <- 0.10
MIN_ABS_RHO  <- 0.50

# search directories for 16S level tables
SEARCH_DIRS <- c(".", "allsample_bundle_from_featuretable")

TAX_LEVELS <- c("phylum", "class", "order", "family", "genus", "species")

# KEEP ONLY THESE GROUPS (drop D*)
KEEP_GROUPS <- c("GR", "MR", "MR+Phage cocktail")

OUTDIR <- "Corr_3groups_noAFMT_16Ssig_q0.1_vs_Metsig_q0.2"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

## ---------- diagnostics + robust paths ----------
cat("[DEBUG] getwd() = ", getwd(), "\n", sep = "")
cat("[DEBUG] Q2.xlsx exists? ", file.exists(file.path(getwd(), MET_QUANT_FILE)), "\n", sep = "")
cat("[DEBUG] stats exists? ", file.exists(file.path(getwd(), MET_STATS_FILE)), "\n", sep = "")

MET_QUANT_XLSX <- normalizePath(MET_QUANT_FILE, winslash = "/", mustWork = TRUE)
MET_STATS_XLSX <- normalizePath(MET_STATS_FILE, winslash = "/", mustWork = TRUE)

## ---------- helpers ----------
infer_group <- function(sample_ids) {
  prefix <- substring(sample_ids, 1, 1)
  grp <- dplyr::case_when(
    prefix == "A" ~ "GR",
    prefix == "B" ~ "MR",
    prefix == "C" ~ "MR+Phage cocktail",
    prefix == "D" ~ "MR+AFMT",
    TRUE ~ NA_character_
  )
  if (any(is.na(grp))) {
    bad <- unique(sample_ids[is.na(grp)])
    stop("Cannot infer group for samples: ", paste(bad, collapse = ", "))
  }
  factor(grp, levels = c("GR","MR","MR+Phage cocktail","MR+AFMT"))
}

keep_samples_3groups <- function(sample_ids) {
  # Keep A*, B*, C* only
  sample_ids[grepl("^[ABC][1-5]$", sample_ids)]
}

read_tsv_loose <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
             stringsAsFactors = FALSE, quote = "", comment.char = "",
             fileEncoding = "UTF-8-BOM")
}

clr_log2 <- function(mat) {
  x <- as.matrix(mat)
  x[is.na(x)] <- 0
  x[x < 0] <- 0
  for (j in seq_len(ncol(x))) {
    v <- x[, j]
    pos <- v[v > 0]
    pc <- if (length(pos) > 0) max(min(pos) / 2, 1e-12) else 1e-6
    v2 <- v + pc
    gm <- exp(mean(log(v2)))
    x[, j] <- log2(v2 / gm)
  }
  x
}

log2_pseudo <- function(mat) {
  x <- as.matrix(mat)
  x[is.na(x)] <- 0
  pos <- x[x > 0]
  pc <- if (length(pos) > 0) max(min(pos) / 2, 1e-12) else 1e-6
  log2(x + pc)
}

safe_wilcox <- function(x, y) {
  out <- tryCatch(stats::wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) 1)
  if (is.na(out) || !is.finite(out)) out <- 1
  out
}

safe_cor <- function(x, y, method = "spearman") {
  out <- tryCatch({
    ct <- suppressWarnings(stats::cor.test(x, y, method = method, exact = FALSE))
    c(rho = unname(ct$estimate), p = ct$p.value)
  }, error = function(e) c(rho = NA_real_, p = 1))
  out
}

find_tax_level_file <- function(level, dirs) {
  patterns <- c(
    paste0("allsample\\.", level, "\\.reabundance\\.xls$"),
    paste0("allsample\\.", level, "\\.abundance\\.xls$"),
    paste0("allsample_", level, "\\.reabundance\\.xls$"),
    paste0("allsample_", level, "\\.abundance\\.xls$")
  )
  hits <- c()
  for (d in dirs) {
    if (!dir.exists(d) && d != ".") next
    for (pat in patterns) {
      cand <- list.files(d, pattern = pat, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
      hits <- c(hits, cand)
    }
  }
  if (length(hits) == 0) return(NA_character_)
  hits <- unique(hits)
  reab <- hits[grepl("reabundance", hits, ignore.case = TRUE)]
  if (length(reab) > 0) return(reab[1])
  hits[1]
}

read_tax_table <- function(path) {
  df <- read_tsv_loose(path) %>% as.data.frame(check.names = FALSE)
  sample_cols <- grep("^[ABCD][1-5]$", names(df), value = TRUE)
  if (length(sample_cols) != 20) {
    stop("16S file ", basename(path), ": cannot detect 20 sample columns A1..D5.")
  }

  # drop D* here
  sample_cols2 <- keep_samples_3groups(sample_cols)
  if (length(sample_cols2) != 15) {
    stop("After dropping D*, expected 15 samples (A/B/C). Detected: ", paste(sample_cols2, collapse=", "))
  }

  tax_col <- names(df)[1]
  taxa <- df[[tax_col]] %>% as.character() %>% str_trim()
  keep <- !is.na(taxa) & taxa != ""
  df <- df[keep, , drop = FALSE]
  taxa <- taxa[keep]

  mat <- df[, sample_cols2, drop = FALSE] %>% as.matrix()
  mat <- apply(mat, 2, function(x) suppressWarnings(as.numeric(x)))
  rownames(mat) <- taxa

  # collapse duplicates by sum
  mat <- rowsum(mat, group = rownames(mat), reorder = FALSE)

  list(mat = mat, samples = sample_cols2)
}

compute_16s_fingerprint <- function(mat, samples, q_cut = 0.1) {
  grp <- infer_group(samples)
  grp <- droplevels(grp[grp %in% KEEP_GROUPS])

  idx_gr <- grp == "GR"
  idx_mr <- grp == "MR"

  clr <- clr_log2(mat)

  pvals <- apply(clr, 1, function(v) safe_wilcox(v[idx_gr], v[idx_mr]))
  qvals <- p.adjust(pvals, method = "BH")

  effect_clr <- rowMeans(clr[, idx_gr, drop = FALSE], na.rm = TRUE) -
    rowMeans(clr[, idx_mr, drop = FALSE], na.rm = TRUE)

  mean_gr <- rowMeans(mat[, idx_gr, drop = FALSE], na.rm = TRUE)
  mean_mr <- rowMeans(mat[, idx_mr, drop = FALSE], na.rm = TRUE)
  pos <- mat[mat > 0]
  pc <- if (length(pos) > 0) max(min(pos) / 2, 1e-12) else 1e-6
  log2fc <- log2((mean_gr + pc) / (mean_mr + pc))

  res <- data.frame(
    feature = rownames(mat),
    effect_CLR_GR_minus_MR = effect_clr,
    log2FC_GR_vs_MR = log2fc,
    p_value = pvals,
    q_value = qvals,
    direction = ifelse(effect_clr >= 0, "GR higher", "MR higher"),
    stringsAsFactors = FALSE
  ) %>% arrange(q_value, desc(abs(effect_CLR_GR_minus_MR)))

  sig <- res %>% filter(!is.na(q_value) & q_value < q_cut)
  list(all = res, sig = sig, clr = clr)
}

run_corr <- function(tax_mat, met_mat, method = "spearman") {
  shared <- intersect(colnames(tax_mat), colnames(met_mat))
  if (length(shared) < 8) stop("Too few shared samples for correlation.")
  X <- tax_mat[, shared, drop = FALSE]
  Y <- met_mat[, shared, drop = FALSE]

  out <- vector("list", nrow(X))
  for (i in seq_len(nrow(X))) {
    x <- as.numeric(X[i, ])
    tmp <- lapply(seq_len(nrow(Y)), function(j) {
      y <- as.numeric(Y[j, ])
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) < 6) return(NULL)
      ct <- safe_cor(x[ok], y[ok], method = method)
      data.frame(
        taxon = rownames(X)[i],
        metabolite = rownames(Y)[j],
        rho = unname(ct["rho"]),
        p_value = unname(ct["p"]),
        n = sum(ok),
        stringsAsFactors = FALSE
      )
    })
    out[[i]] <- bind_rows(tmp)
  }
  res <- bind_rows(out)
  res$q_value <- p.adjust(res$p_value, method = "BH")
  res %>% arrange(q_value, desc(abs(rho)))
}

residualize_by_group <- function(mat, group_factor) {
  X <- as.matrix(mat)
  t(apply(X, 1, function(v) {
    ok <- is.finite(v)
    if (sum(ok) < 6) return(rep(NA_real_, length(v)))
    as.numeric(residuals(lm(v ~ group_factor)))
  }))
}

## ---------- 1) Metabolite fingerprint list (q<0.2, GR vs MR) ----------
stats_raw <- read_excel(MET_STATS_XLSX, sheet = MET_STATS_SHEET) %>%
  as.data.frame(check.names = FALSE)

if (!("metabolite" %in% names(stats_raw))) {
  met_col <- names(stats_raw)[tolower(names(stats_raw)) == "metabolite"][1]
  if (is.na(met_col)) stop("Cannot find 'metabolite' column in stats.")
  names(stats_raw)[names(stats_raw) == met_col] <- "metabolite"
}
stats <- stats_raw
names(stats) <- make.names(names(stats), unique = TRUE)

cand_q <- c("fdr_a_vs_b", "fdr_a_v__b", "fdr_a_v_b")
q_col <- cand_q[cand_q %in% names(stats)][1]
if (is.na(q_col) || q_col == "") q_col <- grep("^fdr.*a.*b$", names(stats), ignore.case = TRUE, value = TRUE)[1]
if (is.na(q_col) || q_col == "") stop("Cannot infer GR vs MR FDR column in stats.")

met_sig <- stats %>%
  transmute(
    metabolite = as.character(metabolite) %>% str_trim(),
    q = suppressWarnings(as.numeric(.data[[q_col]]))
  ) %>%
  filter(!is.na(metabolite) & metabolite != "" & !is.na(q) & q < MET_Q_CUT) %>%
  distinct(metabolite) %>%
  pull(metabolite)

write_csv(data.frame(metabolite = met_sig),
          file.path(OUTDIR, sprintf("MetFingerprint_GRvsMR_q%.2f.csv", MET_Q_CUT)))
cat("[INFO] Met fingerprint q<", MET_Q_CUT, ": n=", length(met_sig), "\n", sep = "")

## ---------- 2) Read metabolomics quant matrix; then DROP D* columns ----------
met_raw <- read_excel(MET_QUANT_XLSX, sheet = excel_sheets(MET_QUANT_XLSX)[1]) %>%
  as.data.frame(check.names = FALSE)

sample_cols_met_all <- grep("^[ABCD][1-5]$", names(met_raw), value = TRUE)
if (length(sample_cols_met_all) != 20) stop("Met quant: cannot detect A1..D5 columns.")

sample_cols_met <- keep_samples_3groups(sample_cols_met_all)
if (length(sample_cols_met) != 15) stop("Met quant: after dropping D*, expected 15 samples A/B/C.")

compound_col <- grep("Compound\\s*Name.*nmol/g", names(met_raw), ignore.case = TRUE, value = TRUE)[1]
if (is.na(compound_col) || compound_col == "") compound_col <- grep("Compound\\s*Name", names(met_raw), ignore.case = TRUE, value = TRUE)[1]
if (is.na(compound_col) || compound_col == "") stop("Met quant: cannot find 'Compound Name' column.")

met_name <- met_raw[[compound_col]] %>% as.character() %>% str_trim()
keep <- !is.na(met_name) & met_name != ""
met_raw2 <- met_raw[keep, , drop = FALSE]
met_name <- met_name[keep]

met_mat <- met_raw2[, sample_cols_met, drop = FALSE] %>% as.matrix()
met_mat <- apply(met_mat, 2, function(x) suppressWarnings(as.numeric(x)))
rownames(met_mat) <- met_name

# collapse duplicates by mean
tmp <- as.data.frame(met_mat) %>% mutate(metabolite = rownames(met_mat))
tmp2 <- tmp %>% group_by(metabolite) %>%
  summarise(across(all_of(sample_cols_met), ~mean(.x, na.rm = TRUE)), .groups="drop")
met_mat <- as.matrix(tmp2[, sample_cols_met, drop = FALSE])
rownames(met_mat) <- tmp2$metabolite

met_sig_in_quant <- intersect(met_sig, rownames(met_mat))
cat("[INFO] Met fingerprint matched in quant (3 groups): ", length(met_sig_in_quant), "\n", sep = "")
met_log2 <- log2_pseudo(met_mat[met_sig_in_quant, , drop = FALSE])
colnames(met_log2) <- sample_cols_met

## ---------- 3) Per 16S level: fingerprint (q<0.1) + correlation (overall + adjusted) ----------
summaries <- list()

for (lvl in TAX_LEVELS) {
  cat("\n[LEVEL] ", lvl, "\n", sep = "")

  f_tax <- find_tax_level_file(lvl, SEARCH_DIRS)
  if (is.na(f_tax)) {
    cat("[SKIP] Cannot find 16S table for level: ", lvl, "\n", sep = "")
    next
  }
  cat("[INFO] Using 16S file: ", f_tax, "\n", sep = "")

  obj <- read_tax_table(f_tax)
  tax_mat <- obj$mat
  samples <- obj$samples  # already A/B/C only
  grp <- infer_group(samples)
  grp <- droplevels(grp[grp %in% KEEP_GROUPS])

  fp <- compute_16s_fingerprint(tax_mat, samples, q_cut = TAX_Q_CUT)

  # save 16S fingerprint tables
  out_all <- file.path(OUTDIR, sprintf("16S_%s_GRvsMR_ALL_3groups_noAFMT.csv", lvl))
  out_sig <- file.path(OUTDIR, sprintf("16S_%s_GRvsMR_Fingerprint_q%.2f_3groups_noAFMT.csv", lvl, TAX_Q_CUT))
  write_csv(fp$all, out_all)
  write_csv(fp$sig, out_sig)

  cat("[INFO] 16S fingerprint q<", TAX_Q_CUT, ": n=", nrow(fp$sig), "\n", sep = "")

  if (nrow(fp$sig) == 0 || length(met_sig_in_quant) == 0) {
    cat("[SKIP] No fingerprint taxa or no metabolite fingerprint matched.\n")
    next
  }

  sig_taxa <- fp$sig$feature
  tax_clr_sig <- fp$clr[sig_taxa, , drop = FALSE]
  colnames(tax_clr_sig) <- samples

  shared <- intersect(colnames(tax_clr_sig), colnames(met_log2))
  if (length(shared) < 8) stop("Too few shared samples at level ", lvl)

  # overall
  corr_overall <- run_corr(tax_clr_sig, met_log2, method = METHOD_COR) %>%
    mutate(level = lvl, mode = "overall")

  # group-adjusted
  tax_resid <- residualize_by_group(tax_clr_sig[, shared, drop = FALSE], grp[match(shared, samples)])
  rownames(tax_resid) <- rownames(tax_clr_sig); colnames(tax_resid) <- shared

  met_resid <- residualize_by_group(met_log2[, shared, drop = FALSE], grp[match(shared, samples)])
  rownames(met_resid) <- rownames(met_log2); colnames(met_resid) <- shared

  corr_adj <- run_corr(tax_resid, met_resid, method = METHOD_COR) %>%
    mutate(level = lvl, mode = "group_adjusted")

  corr_all <- bind_rows(corr_overall, corr_adj) %>% arrange(mode, q_value, desc(abs(rho)))
  out_corr <- file.path(OUTDIR, sprintf("CORR_%s_16Sfp_q%.2f__Metfp_q%.2f_Spearman_3groups_noAFMT.csv",
                                        lvl, TAX_Q_CUT, MET_Q_CUT))
  write_csv(corr_all, out_corr)

  edges <- corr_all %>%
    filter(!is.na(rho) & !is.na(q_value) & q_value < CORR_Q_CUT & abs(rho) >= MIN_ABS_RHO) %>%
    arrange(mode, q_value, desc(abs(rho)))
  out_edges <- file.path(OUTDIR, sprintf("EDGES_%s_q%.2f_rho%.2f_3groups_noAFMT.csv",
                                         lvl, CORR_Q_CUT, MIN_ABS_RHO))
  write_csv(edges, out_edges)

  cat("[OK] Corr saved: ", out_corr, "\n", sep = "")
  cat("[OK] Filtered edges: n=", nrow(edges), " -> ", out_edges, "\n", sep = "")

  if (HAS_PHEATMAP && nrow(edges) > 1) {
    for (md in unique(edges$mode)) {
      e2 <- edges %>% filter(mode == md)
      if (nrow(e2) < 2) next
      mat_rho <- e2 %>%
        select(taxon, metabolite, rho) %>%
        distinct() %>%
        pivot_wider(names_from = metabolite, values_from = rho, values_fill = 0) %>%
        as.data.frame()
      rownames(mat_rho) <- mat_rho$taxon
      mat_rho$taxon <- NULL
      mat_rho <- as.matrix(mat_rho)

      pdf(file.path(OUTDIR, sprintf("HEATMAP_%s_%s_q%.2f_rho%.2f_3groups_noAFMT.pdf",
                                    lvl, md, CORR_Q_CUT, MIN_ABS_RHO)),
          width = 10, height = max(5, 0.22 * nrow(mat_rho) + 2))
      pheatmap::pheatmap(mat_rho,
                         cluster_rows = TRUE,
                         cluster_cols = TRUE,
                         main = "",
                         fontsize_row = 7,
                         fontsize_col = 7,
                         border_color = NA)
      dev.off()
    }
  }

  summaries[[lvl]] <- data.frame(
    level = lvl,
    n_tax_fp = nrow(fp$sig),
    n_met_fp = length(met_sig_in_quant),
    n_edges_overall = sum(edges$mode == "overall"),
    n_edges_adj = sum(edges$mode == "group_adjusted"),
    stringsAsFactors = FALSE
  )
}

if (length(summaries) > 0) {
  summ <- bind_rows(summaries)
  write_csv(summ, file.path(OUTDIR, "SUMMARY_levels_3groups_noAFMT.csv"))
  cat("\n[SUMMARY]\n")
  print(summ)
}

cat("\n[DONE] Outputs in: ", OUTDIR, "\n", sep = "")
