#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(vegan)
})

# ============================================================
# 3-group metabolomics: PCA + PERMANOVA/betadisper + Dist2GR
# Sankey/Venn-aligned rewrite
#
# Groups kept:
#   A -> GR
#   B -> MR
#   C -> MR+Phage cocktail
#
# Inputs:
#   - Q2.xlsx (sheet: sample)
#   - Sample_Metadata_3groups.csv
#   - 3-group stats table from the rewritten stats script
#
# Main outputs:
#   results_3group_metabolomics_PCA_PERMANOVA_Dist2GR_sankey_aligned/
#     01_input_checks/
#     02_PCA_plots/
#     03_PCA_tables/
#     04_PERMANOVA_betadisper/
#     05_Dist2GR/
#
# Processing:
#   - per-metabolite 0 / NA -> half of minimum positive value
#   - log2 transform
#   - row z-score
#   - PCA on samples
#   - Euclidean distance in transformed feature space
#   - PERMANOVA + betadisper on 3 groups
#   - AB signature defined by fdr_A_vs_B < 0.2 from current 3-group stats table
# ============================================================

options(stringsAsFactors = FALSE)

# ----------------------------
# User parameters
# ----------------------------
MET_QUANT_XLSX  <- "Q2.xlsx"
MET_QUANT_SHEET <- "sample"

META_CSV <- "Sample_Metadata_3groups.csv"

# Prefer the 3-group aligned stats generated previously
STATS_XLSX <- file.path(
  "results_3group_metabolomics_stats_volcano_heatmap_sankey_aligned",
  "01_stats",
  "metabolomics_stats_all_results_3groups_sankey_aligned.xlsx"
)
STATS_CSV_FALLBACK <- file.path(
  "results_3group_metabolomics_stats_volcano_heatmap_sankey_aligned",
  "01_stats",
  "metabolomics_stats_all_results_3groups_sankey_aligned.csv"
)
STATS_SHEET <- "all_results"

OUTDIR <- "results_3group_metabolomics_PCA_PERMANOVA_Dist2GR_sankey_aligned"

GROUP_LEVELS <- c("GR", "MR", "MR+Phage cocktail")
GROUP_COLORS <- c(
  "GR" = "#7CAE00",
  "MR" = "#F8766D",
  "MR+Phage cocktail" = "#00BFC4"
)

FINGERPRINT_Q <- 0.20
TOP_N_FALLBACK <- 20
PC_N <- 5
ELLIPSE_LEVEL <- 0.68
LABEL_SAMPLES <- FALSE
N_PERM <- 9999

# Dist2GR focus comparisons to annotate on plots
FOCUS_PAIRS <- list(
  c("MR", "MR+Phage cocktail")
)

# ----------------------------
# Directories
# ----------------------------
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "01_input_checks"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "02_PCA_plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "03_PCA_tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "04_PERMANOVA_betadisper"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "05_Dist2GR"), showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Helpers
# ----------------------------
msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n")
}

stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
}

safe_write_csv <- function(x, file) {
  readr::write_csv(x, file, na = "")
}

pick_first_col <- function(nms, candidates) {
  hit <- intersect(candidates, nms)
  if (length(hit) > 0) return(hit[1])
  NA_character_
}

pair_label <- function(g1, g2) paste0(g1, "_vs_", g2)

impute_row_halfmin <- function(v) {
  v <- suppressWarnings(as.numeric(v))
  if (all(is.na(v))) return(rep(NA_real_, length(v)))
  min_pos <- suppressWarnings(min(v[v > 0], na.rm = TRUE))
  if (!is.finite(min_pos)) return(rep(NA_real_, length(v)))
  v[is.na(v) | v == 0] <- min_pos / 2
  v
}

row_zscore <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sdv <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) return(rep(0, length(x)))
  (x - mu) / sdv
}

read_stats_table <- function(xlsx_path, csv_fallback, sheet) {
  if (file.exists(xlsx_path)) {
    msg("Reading stats xlsx: ", xlsx_path)
    return(readxl::read_excel(xlsx_path, sheet = sheet))
  }
  if (file.exists(csv_fallback)) {
    msg("Reading stats csv fallback: ", csv_fallback)
    return(readr::read_csv(csv_fallback, show_col_types = FALSE))
  }
  stop("Cannot find stats table:\n - ", xlsx_path, "\n - ", csv_fallback)
}

pairwise_wilcox_3groups <- function(df, value_col, group_col = "group") {
  df <- df %>%
    filter(!is.na(.data[[value_col]])) %>%
    mutate(group = factor(.data[[group_col]], levels = GROUP_LEVELS))

  pairs <- combn(GROUP_LEVELS, 2, simplify = FALSE)

  out <- bind_rows(lapply(pairs, function(p) {
    g1 <- p[1]; g2 <- p[2]
    x <- df %>% filter(group == g1) %>% pull(.data[[value_col]])
    y <- df %>% filter(group == g2) %>% pull(.data[[value_col]])
    wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
    tibble(
      group1 = g1,
      group2 = g2,
      n1 = length(x),
      n2 = length(y),
      median1 = median(x, na.rm = TRUE),
      median2 = median(y, na.rm = TRUE),
      delta_median = median1 - median2,
      stat_W = unname(wt$statistic),
      p_value = wt$p.value
    )
  })) %>%
    mutate(q_value = p.adjust(p_value, method = "BH"))
  out
}

dist_to_GR_centroid_scores <- function(scores_df, pcs) {
  stopifnot(all(pcs %in% names(scores_df)))
  cent <- scores_df %>%
    filter(group == "GR") %>%
    summarise(across(all_of(pcs), mean, na.rm = TRUE))
  center <- as.numeric(cent[1, pcs, drop = TRUE])

  mat <- as.matrix(scores_df[, pcs, drop = FALSE])
  d <- sqrt(rowSums((mat - matrix(center, nrow = nrow(mat), ncol = length(center), byrow = TRUE))^2))

  tibble(
    sample = scores_df$sample,
    group = scores_df$group,
    pcs_used = paste(pcs, collapse = "+"),
    dist_to_GR = as.numeric(d)
  )
}

dist_to_GR_centroid_feature_space <- function(Y, samples, groups) {
  # Y = samples x features matrix
  gr_idx <- which(groups == "GR")
  if (length(gr_idx) < 2) stop("Need at least 2 GR samples to compute centroid.")
  gr_cent <- colMeans(Y[gr_idx, , drop = FALSE], na.rm = TRUE)
  d <- sqrt(rowSums((Y - matrix(gr_cent, nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE))^2))
  tibble(
    sample = samples,
    group = factor(groups, levels = GROUP_LEVELS),
    dist_to_GR = as.numeric(d)
  )
}

make_dist2gr_plot <- function(dist_df, focus_pairs, ylab) {
  ymax <- max(dist_df$dist_to_GR, na.rm = TRUE)
  ymin <- min(dist_df$dist_to_GR, na.rm = TRUE)
  step <- (ymax - ymin) * 0.12
  if (!is.finite(step) || step == 0) step <- 0.5

  x_pos <- setNames(seq_along(GROUP_LEVELS), GROUP_LEVELS)

  focus_tests <- bind_rows(lapply(focus_pairs, function(p) {
    g1 <- p[1]; g2 <- p[2]
    x <- dist_df %>% filter(group == g1) %>% pull(dist_to_GR)
    y <- dist_df %>% filter(group == g2) %>% pull(dist_to_GR)
    wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
    tibble(
      group1 = g1,
      group2 = g2,
      p_value = wt$p.value
    )
  })) %>%
    mutate(q_value = p.adjust(p_value, method = "BH")) %>%
    mutate(
      label = paste0("p=", signif(p_value, 3), "\nq=", signif(q_value, 3)),
      x1 = x_pos[group1],
      x2 = x_pos[group2],
      y = ymax + row_number() * step
    )

  p <- ggplot(dist_df, aes(x = group, y = dist_to_GR, color = group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.7) +
    geom_jitter(width = 0.12, size = 2.2, alpha = 0.85) +
    scale_color_manual(values = GROUP_COLORS) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12),
      axis.title.y = element_text(size = 13)
    ) +
    labs(x = NULL, y = ylab) +
    coord_cartesian(ylim = c(ymin, ymax + (nrow(focus_tests) + 2) * step))

  if (nrow(focus_tests) > 0) {
    p <- p +
      geom_segment(data = focus_tests,
                   aes(x = x1, xend = x2, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = 0.7, color = "black") +
      geom_segment(data = focus_tests,
                   aes(x = x1, xend = x1, y = y, yend = y - step * 0.15),
                   inherit.aes = FALSE, linewidth = 0.7, color = "black") +
      geom_segment(data = focus_tests,
                   aes(x = x2, xend = x2, y = y, yend = y - step * 0.15),
                   inherit.aes = FALSE, linewidth = 0.7, color = "black") +
      geom_text(data = focus_tests,
                aes(x = (x1 + x2) / 2, y = y + step * 0.05, label = label),
                inherit.aes = FALSE, size = 3.6, color = "black")
  }

  list(plot = p, focus_tests = focus_tests)
}

preprocess_quant_matrix <- function(raw_df, sample_cols, feature_keep = NULL, feature_name_col = "metabolite") {
  df <- raw_df %>%
    mutate(metabolite = as.character(.data[[feature_name_col]])) %>%
    select(metabolite, all_of(sample_cols))

  if (!is.null(feature_keep)) {
    df <- df %>% filter(metabolite %in% feature_keep)
  }

  df2 <- df %>%
    rowwise() %>%
    mutate(across(all_of(sample_cols), ~ impute_row_halfmin(.x))) %>%
    ungroup() %>%
    mutate(across(all_of(sample_cols), ~ log2(.x)))

  mat <- df2 %>% as.data.frame()
  rownames(mat) <- make.unique(mat$metabolite)
  mat$metabolite <- NULL
  X <- as.matrix(mat)

  keep <- apply(X, 1, function(x) all(is.finite(x)) && sd(x) > 0)
  X <- X[keep, , drop = FALSE]

  Xz <- t(apply(X, 1, row_zscore))
  Xz[!is.finite(Xz)] <- 0

  list(
    log2_matrix = X,
    z_matrix = Xz,
    kept_features = rownames(Xz)
  )
}

run_pca_from_quant <- function(raw_df, sample_cols, meta_df, feature_keep = NULL, feature_name_col = "metabolite", pc_n = 5) {
  prep <- preprocess_quant_matrix(raw_df, sample_cols, feature_keep = feature_keep, feature_name_col = feature_name_col)
  Xz <- prep$z_matrix

  pca <- prcomp(t(Xz), center = FALSE, scale. = FALSE)

  k <- min(pc_n, ncol(pca$x))
  scores <- as.data.frame(pca$x[, 1:k, drop = FALSE]) %>%
    rownames_to_column("sample") %>%
    left_join(meta_df, by = "sample") %>%
    rename(group = group)

  eig <- pca$sdev^2
  variance <- tibble(
    PC = paste0("PC", seq_along(eig)),
    sdev = as.numeric(pca$sdev),
    eigenvalue = as.numeric(eig),
    variance_pct = as.numeric(eig / sum(eig) * 100)
  )

  loadings <- as.data.frame(pca$rotation[, 1:k, drop = FALSE]) %>%
    rownames_to_column("metabolite")

  list(
    prep = prep,
    pca = pca,
    scores = scores,
    variance = variance,
    loadings = loadings
  )
}

make_pca_plot <- function(scores_df, variance_df, label_samples = FALSE, ellipse_level = 0.68) {
  pc1_lab <- paste0("PC1 (", round(variance_df$variance_pct[1], 1), "%)")
  pc2_lab <- paste0("PC2 (", round(variance_df$variance_pct[2], 1), "%)")

  p <- ggplot(scores_df, aes(PC1, PC2, color = group, fill = group)) +
    stat_ellipse(geom = "polygon", type = "t", level = ellipse_level, alpha = 0.18, color = NA) +
    stat_ellipse(type = "t", level = ellipse_level, linewidth = 0.6) +
    geom_point(size = 3) +
    scale_color_manual(values = GROUP_COLORS) +
    scale_fill_manual(values = GROUP_COLORS) +
    labs(x = pc1_lab, y = pc2_lab) +
    theme_classic(base_size = 14) +
    theme(legend.title = element_blank())

  if (label_samples) {
    p <- p + geom_text(aes(label = sample), size = 3, vjust = -0.6)
  }
  p
}

run_feature_space_perm <- function(Xz, meta_df, tag, outdir, n_perm = 9999) {
  Y <- t(Xz)  # samples x features
  meta <- meta_df %>%
    filter(sample %in% rownames(Y)) %>%
    mutate(group = factor(group, levels = GROUP_LEVELS))
  Y <- Y[meta$sample, , drop = FALSE]

  d <- dist(Y, method = "euclidean")

  ad <- adonis2(d ~ group, data = meta, permutations = n_perm)
  ad_df <- as.data.frame(ad)
  ad_df$term <- rownames(ad_df)
  rownames(ad_df) <- NULL
  ad_df <- ad_df %>% select(term, everything())
  safe_write_csv(ad_df, file.path(outdir, paste0("PERMANOVA_global_", tag, ".csv")))

  pairs <- combn(GROUP_LEVELS, 2, simplify = FALSE)
  pair_res <- bind_rows(lapply(pairs, function(p) {
    g1 <- p[1]; g2 <- p[2]
    idx <- meta$group %in% c(g1, g2)
    meta2 <- droplevels(meta[idx, , drop = FALSE])
    Y2 <- Y[idx, , drop = FALSE]
    d2 <- dist(Y2, method = "euclidean")
    ad2 <- adonis2(d2 ~ group, data = meta2, permutations = n_perm)
    a <- as.data.frame(ad2)
    tibble(
      group1 = g1,
      group2 = g2,
      R2 = a$R2[1],
      F = a$F[1],
      p_value = a$`Pr(>F)`[1]
    )
  })) %>%
    mutate(q_value = p.adjust(p_value, method = "BH")) %>%
    arrange(p_value)
  safe_write_csv(pair_res, file.path(outdir, paste0("PERMANOVA_pairwise_", tag, ".csv")))

  bd <- betadisper(d, meta$group)
  bd_aov <- anova(bd)
  bd_perm <- permutest(bd, permutations = n_perm)

  capture.output(
    list(
      "ANOVA(betadisper):", bd_aov, "",
      "permutest(betadisper):", bd_perm, "",
      "TukeyHSD(betadisper):", TukeyHSD(bd)
    ),
    file = file.path(outdir, paste0("betadisper_", tag, ".txt"))
  )

  disp_df <- data.frame(
    sample = meta$sample,
    group = meta$group,
    dist_to_centroid = bd$distances
  )
  safe_write_csv(disp_df, file.path(outdir, paste0("betadisper_distances_", tag, ".csv")))

  p_disp <- ggplot(disp_df, aes(x = group, y = dist_to_centroid, color = group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.7) +
    geom_jitter(width = 0.12, size = 2.2, alpha = 0.85) +
    scale_color_manual(values = GROUP_COLORS) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") +
    labs(x = NULL, y = "Distance to group centroid (dispersion)")

  ggsave(file.path(outdir, paste0("betadisper_boxplot_", tag, ".png")), p_disp, width = 7.6, height = 5.0, dpi = 300)
  ggsave(file.path(outdir, paste0("betadisper_boxplot_", tag, ".pdf")), p_disp, width = 7.6, height = 5.0)

  pcoa <- cmdscale(d, k = 2, eig = TRUE)
  pts <- as.data.frame(pcoa$points)
  colnames(pts) <- c("Axis1", "Axis2")
  pts$sample <- rownames(pts)
  pts <- pts %>% left_join(meta, by = "sample")

  var1 <- round(100 * pcoa$eig[1] / sum(pcoa$eig), 1)
  var2 <- round(100 * pcoa$eig[2] / sum(pcoa$eig), 1)

  p_pcoa <- ggplot(pts, aes(Axis1, Axis2, color = group, fill = group)) +
    stat_ellipse(geom = "polygon", type = "t", level = ELLIPSE_LEVEL, alpha = 0.18, color = NA) +
    stat_ellipse(type = "t", level = ELLIPSE_LEVEL, linewidth = 0.6) +
    geom_point(size = 3) +
    scale_color_manual(values = GROUP_COLORS) +
    scale_fill_manual(values = GROUP_COLORS) +
    theme_classic(base_size = 14) +
    theme(legend.title = element_blank()) +
    labs(x = paste0("PCoA1 (", var1, "%)"), y = paste0("PCoA2 (", var2, "%)"))

  ggsave(file.path(outdir, paste0("PCoA_", tag, ".png")), p_pcoa, width = 7.2, height = 5.2, dpi = 300)
  ggsave(file.path(outdir, paste0("PCoA_", tag, ".pdf")), p_pcoa, width = 7.2, height = 5.2)

  invisible(list(
    permanova_global = ad_df,
    permanova_pairwise = pair_res,
    betadisper_distances = disp_df
  ))
}

# ----------------------------
# 1) Read metadata
# ----------------------------
msg("Reading metadata...")
stop_if_missing(META_CSV)
meta <- readr::read_csv(META_CSV, show_col_types = FALSE)

if (!all(c("sample", "group") %in% names(meta))) {
  stop("Sample_Metadata_3groups.csv must contain columns: sample, group")
}

meta <- meta %>%
  filter(group %in% GROUP_LEVELS) %>%
  mutate(group = factor(group, levels = GROUP_LEVELS)) %>%
  arrange(group, sample)

sample_order <- meta$sample
safe_write_csv(meta, file.path(OUTDIR, "01_input_checks", "metadata_used.csv"))

# ----------------------------
# 2) Read quantitative table
# ----------------------------
msg("Reading quantitative metabolomics table...")
stop_if_missing(MET_QUANT_XLSX)

raw <- readxl::read_excel(MET_QUANT_XLSX, sheet = MET_QUANT_SHEET)
raw <- as.data.frame(raw)

sample_cols <- intersect(names(raw), sample_order)
if (length(sample_cols) == 0) {
  sample_cols <- names(raw)[str_detect(names(raw), "^[ABC]\\d+$")]
}
if (length(sample_cols) == 0) stop("Cannot detect 3-group sample columns in Q2.xlsx.")

# enforce metadata order
sample_cols <- sample_order[sample_order %in% sample_cols]

met_col_raw <- pick_first_col(
  names(raw),
  c("metabolite", "Compound Name (nmol/g)", "Compound Name", "代谢物", "Metabolite", "Name")
)
if (is.na(met_col_raw)) met_col_raw <- names(raw)[1]

raw2 <- raw %>%
  mutate(metabolite = as.character(.data[[met_col_raw]])) %>%
  select(metabolite, all_of(sample_cols))

safe_write_csv(
  tibble(sample = sample_cols, group = meta$group[match(sample_cols, meta$sample)]),
  file.path(OUTDIR, "01_input_checks", "sample_columns_used.csv")
)

safe_write_csv(
  tibble(metabolite = raw2$metabolite),
  file.path(OUTDIR, "01_input_checks", "metabolite_names_in_quant_table.csv")
)

# ----------------------------
# 3) Read 3-group stats and define AB signature
# ----------------------------
msg("Reading 3-group stats table...")
stats <- read_stats_table(STATS_XLSX, STATS_CSV_FALLBACK, STATS_SHEET)
stats <- as.data.frame(stats)

need_stats_cols <- c("metabolite", "fdr_A_vs_B")
miss <- setdiff(need_stats_cols, names(stats))
if (length(miss) > 0) stop("Stats table missing required columns: ", paste(miss, collapse = ", "))

fingerprint_mets <- stats %>%
  transmute(
    metabolite = as.character(metabolite),
    fdr_A_vs_B = suppressWarnings(as.numeric(fdr_A_vs_B))
  ) %>%
  filter(!is.na(metabolite), metabolite != "", !is.na(fdr_A_vs_B), fdr_A_vs_B < FINGERPRINT_Q) %>%
  arrange(fdr_A_vs_B) %>%
  distinct(metabolite) %>%
  pull(metabolite)

if (length(fingerprint_mets) < 2) {
  fingerprint_mets <- stats %>%
    transmute(
      metabolite = as.character(metabolite),
      fdr_A_vs_B = suppressWarnings(as.numeric(fdr_A_vs_B))
    ) %>%
    filter(!is.na(metabolite), metabolite != "", !is.na(fdr_A_vs_B)) %>%
    arrange(fdr_A_vs_B) %>%
    slice(1:TOP_N_FALLBACK) %>%
    distinct(metabolite) %>%
    pull(metabolite)
  msg("AB signature < 2, fallback to Top ", TOP_N_FALLBACK, " by fdr_A_vs_B.")
} else {
  msg("Using AB signature (q<", FINGERPRINT_Q, "), n = ", length(fingerprint_mets))
}

n_match <- sum(fingerprint_mets %in% raw2$metabolite)
if (n_match < length(fingerprint_mets)) {
  stop(
    "Raw quant table metabolite names do not sufficiently match stats table metabolite names: ",
    n_match, "/", length(fingerprint_mets)
  )
}

safe_write_csv(
  tibble(metabolite = fingerprint_mets),
  file.path(OUTDIR, "01_input_checks", paste0("fingerprint_metabolites_q", FINGERPRINT_Q, ".csv"))
)

# ----------------------------
# 4) PCA: all metabolites and AB signature
# ----------------------------
msg("Running PCA...")

res_all <- run_pca_from_quant(
  raw_df = raw2,
  sample_cols = sample_cols,
  meta_df = meta,
  feature_keep = NULL,
  feature_name_col = "metabolite",
  pc_n = PC_N
)

res_fp <- run_pca_from_quant(
  raw_df = raw2,
  sample_cols = sample_cols,
  meta_df = meta,
  feature_keep = fingerprint_mets,
  feature_name_col = "metabolite",
  pc_n = PC_N
)

# Save transformed matrices
safe_write_csv(
  as.data.frame(res_all$prep$log2_matrix) %>% rownames_to_column("metabolite"),
  file.path(OUTDIR, "03_PCA_tables", "AllMet_log2_matrix.csv")
)
safe_write_csv(
  as.data.frame(res_all$prep$z_matrix) %>% rownames_to_column("metabolite"),
  file.path(OUTDIR, "03_PCA_tables", "AllMet_row_z_matrix.csv")
)
safe_write_csv(
  as.data.frame(res_fp$prep$log2_matrix) %>% rownames_to_column("metabolite"),
  file.path(OUTDIR, "03_PCA_tables", paste0("Fingerprint_q", FINGERPRINT_Q, "_log2_matrix.csv"))
)
safe_write_csv(
  as.data.frame(res_fp$prep$z_matrix) %>% rownames_to_column("metabolite"),
  file.path(OUTDIR, "03_PCA_tables", paste0("Fingerprint_q", FINGERPRINT_Q, "_row_z_matrix.csv"))
)

# PCA plots
p_all <- make_pca_plot(res_all$scores, res_all$variance, label_samples = LABEL_SAMPLES, ellipse_level = ELLIPSE_LEVEL)
p_fp  <- make_pca_plot(res_fp$scores, res_fp$variance, label_samples = LABEL_SAMPLES, ellipse_level = ELLIPSE_LEVEL)

ggsave(file.path(OUTDIR, "02_PCA_plots", "PCA_AllMet_ellipse_3groups.png"), p_all, width = 7.2, height = 5.2, dpi = 300)
ggsave(file.path(OUTDIR, "02_PCA_plots", "PCA_AllMet_ellipse_3groups.pdf"), p_all, width = 7.2, height = 5.2)

ggsave(file.path(OUTDIR, "02_PCA_plots", paste0("PCA_ABsig_q", FINGERPRINT_Q, "_ellipse_3groups.png")), p_fp, width = 7.2, height = 5.2, dpi = 300)
ggsave(file.path(OUTDIR, "02_PCA_plots", paste0("PCA_ABsig_q", FINGERPRINT_Q, "_ellipse_3groups.pdf")), p_fp, width = 7.2, height = 5.2)

# PCA tables export helper
export_one_pca <- function(res, tag) {
  safe_write_csv(res$scores,   file.path(OUTDIR, "03_PCA_tables", paste0("PCA_", tag, "_scores_PC1toPC", PC_N, ".csv")))
  safe_write_csv(res$variance, file.path(OUTDIR, "03_PCA_tables", paste0("PCA_", tag, "_variance.csv")))
  safe_write_csv(res$loadings, file.path(OUTDIR, "03_PCA_tables", paste0("PCA_", tag, "_loadings_PC1toPC", PC_N, ".csv")))

  pcs <- paste0("PC", 1:min(PC_N, ncol(res$scores) - 2))

  cent <- res$scores %>%
    group_by(group) %>%
    summarise(across(all_of(pcs), mean, na.rm = TRUE), .groups = "drop")
  safe_write_csv(cent, file.path(OUTDIR, "03_PCA_tables", paste0("PCA_", tag, "_group_centroids.csv")))

  pw_pc <- bind_rows(lapply(pcs, function(pc) {
    pairwise_wilcox_3groups(res$scores, value_col = pc, group_col = "group") %>%
      mutate(dataset = tag, PC = pc) %>%
      select(dataset, PC, everything())
  }))
  safe_write_csv(pw_pc, file.path(OUTDIR, "03_PCA_tables", paste0("PCA_", tag, "_pairwise_Wilcoxon_3groups_PC_scores.csv")))

  d12 <- dist_to_GR_centroid_scores(res$scores, pcs = intersect(c("PC1", "PC2"), pcs))
  d15 <- dist_to_GR_centroid_scores(res$scores, pcs = pcs)

  safe_write_csv(d12, file.path(OUTDIR, "03_PCA_tables", paste0("PCA_", tag, "_dist_to_GR_PC1_2_perSample.csv")))
  safe_write_csv(d15, file.path(OUTDIR, "03_PCA_tables", paste0("PCA_", tag, "_dist_to_GR_PC1_5_perSample.csv")))

  pw_d12 <- pairwise_wilcox_3groups(d12, value_col = "dist_to_GR", group_col = "group") %>%
    mutate(dataset = tag, pcs_used = "PC1+PC2") %>%
    select(dataset, pcs_used, everything())
  pw_d15 <- pairwise_wilcox_3groups(d15, value_col = "dist_to_GR", group_col = "group") %>%
    mutate(dataset = tag, pcs_used = paste(pcs, collapse = "+")) %>%
    select(dataset, pcs_used, everything())

  safe_write_csv(pw_d12, file.path(OUTDIR, "03_PCA_tables", paste0("PCA_", tag, "_pairwise_Wilcoxon_3groups_dist_to_GR_PC1_2.csv")))
  safe_write_csv(pw_d15, file.path(OUTDIR, "03_PCA_tables", paste0("PCA_", tag, "_pairwise_Wilcoxon_3groups_dist_to_GR_PC1_5.csv")))
}

export_one_pca(res_all, "AllMet")
export_one_pca(res_fp, paste0("Fingerprint_q", FINGERPRINT_Q))

# ----------------------------
# 5) PERMANOVA + betadisper
# ----------------------------
msg("Running PERMANOVA + betadisper...")

run_feature_space_perm(
  Xz = res_all$prep$z_matrix,
  meta_df = meta,
  tag = "AllMet_featureSpace_3groups",
  outdir = file.path(OUTDIR, "04_PERMANOVA_betadisper"),
  n_perm = N_PERM
)

run_feature_space_perm(
  Xz = res_fp$prep$z_matrix,
  meta_df = meta,
  tag = paste0("ABsig_q", FINGERPRINT_Q, "_featureSpace_3groups"),
  outdir = file.path(OUTDIR, "04_PERMANOVA_betadisper"),
  n_perm = N_PERM
)

# ----------------------------
# 6) Dist2GR plots and tables
# ----------------------------
msg("Running Dist2GR analyses...")

# 6.1 AB signature PCA PC1+PC2
dist_pc12 <- dist_to_GR_centroid_scores(res_fp$scores, pcs = c("PC1", "PC2"))
safe_write_csv(dist_pc12, file.path(OUTDIR, "05_Dist2GR", paste0("Dist2GR_ABsig_q", FINGERPRINT_Q, "_PC1_to_PC2_per_sample.csv")))

summary_pc12 <- dist_pc12 %>%
  group_by(group) %>%
  summarise(
    n = n(),
    median = median(dist_to_GR),
    q1 = quantile(dist_to_GR, 0.25),
    q3 = quantile(dist_to_GR, 0.75),
    .groups = "drop"
  )
safe_write_csv(summary_pc12, file.path(OUTDIR, "05_Dist2GR", paste0("Dist2GR_ABsig_q", FINGERPRINT_Q, "_PC1_to_PC2_group_summary.csv")))

tests_pc12 <- pairwise_wilcox_3groups(dist_pc12, value_col = "dist_to_GR", group_col = "group")
safe_write_csv(tests_pc12, file.path(OUTDIR, "05_Dist2GR", paste0("Dist2GR_ABsig_q", FINGERPRINT_Q, "_PC1_to_PC2_pairwise_tests.csv")))

plot_pc12 <- make_dist2gr_plot(
  dist_df = dist_pc12,
  focus_pairs = FOCUS_PAIRS,
  ylab = "Distance to GR centroid (PC1+PC2)"
)

ggsave(file.path(OUTDIR, "05_Dist2GR", paste0("Dist2GR_ABsig_q", FINGERPRINT_Q, "_PC1_to_PC2.png")),
       plot_pc12$plot, width = 7.6, height = 5.2, dpi = 300)
ggsave(file.path(OUTDIR, "05_Dist2GR", paste0("Dist2GR_ABsig_q", FINGERPRINT_Q, "_PC1_to_PC2.pdf")),
       plot_pc12$plot, width = 7.6, height = 5.2)

# 6.2 AB signature PCA PC1..PC5
pcs_1to5 <- paste0("PC", 1:min(PC_N, ncol(res_fp$scores) - 2))
dist_pc15 <- dist_to_GR_centroid_scores(res_fp$scores, pcs = pcs_1to5)
safe_write_csv(dist_pc15, file.path(OUTDIR, "05_Dist2GR", paste0("Dist2GR_ABsig_q", FINGERPRINT_Q, "_PC1_to_PC5_per_sample.csv")))

summary_pc15 <- dist_pc15 %>%
  group_by(group) %>%
  summarise(
    n = n(),
    median = median(dist_to_GR),
    q1 = quantile(dist_to_GR, 0.25),
    q3 = quantile(dist_to_GR, 0.75),
    .groups = "drop"
  )
safe_write_csv(summary_pc15, file.path(OUTDIR, "05_Dist2GR", paste0("Dist2GR_ABsig_q", FINGERPRINT_Q, "_PC1_to_PC5_group_summary.csv")))

tests_pc15 <- pairwise_wilcox_3groups(dist_pc15, value_col = "dist_to_GR", group_col = "group")
safe_write_csv(tests_pc15, file.path(OUTDIR, "05_Dist2GR", paste0("Dist2GR_ABsig_q", FINGERPRINT_Q, "_PC1_to_PC5_pairwise_tests.csv")))

plot_pc15 <- make_dist2gr_plot(
  dist_df = dist_pc15,
  focus_pairs = FOCUS_PAIRS,
  ylab = paste0("Distance to GR centroid (", paste(pcs_1to5, collapse = "+"), ")")
)

ggsave(file.path(OUTDIR, "05_Dist2GR", paste0("Dist2GR_ABsig_q", FINGERPRINT_Q, "_PC1_to_PC5.png")),
       plot_pc15$plot, width = 7.6, height = 5.2, dpi = 300)
ggsave(file.path(OUTDIR, "05_Dist2GR", paste0("Dist2GR_ABsig_q", FINGERPRINT_Q, "_PC1_to_PC5.pdf")),
       plot_pc15$plot, width = 7.6, height = 5.2)

# 6.3 AB signature feature space direct distance
Y_fp <- t(res_fp$prep$z_matrix)
dist_feature <- dist_to_GR_centroid_feature_space(
  Y = Y_fp,
  samples = rownames(Y_fp),
  groups = as.character(meta$group[match(rownames(Y_fp), meta$sample)])
)

safe_write_csv(dist_feature, file.path(OUTDIR, "05_Dist2GR", paste0("Dist2GR_ABsig_q", FINGERPRINT_Q, "_featureSpace_per_sample.csv")))

summary_feature <- dist_feature %>%
  group_by(group) %>%
  summarise(
    n = n(),
    median = median(dist_to_GR),
    q1 = quantile(dist_to_GR, 0.25),
    q3 = quantile(dist_to_GR, 0.75),
    .groups = "drop"
  )
safe_write_csv(summary_feature, file.path(OUTDIR, "05_Dist2GR", paste0("Dist2GR_ABsig_q", FINGERPRINT_Q, "_featureSpace_group_summary.csv")))

tests_feature <- pairwise_wilcox_3groups(dist_feature, value_col = "dist_to_GR", group_col = "group")
safe_write_csv(tests_feature, file.path(OUTDIR, "05_Dist2GR", paste0("Dist2GR_ABsig_q", FINGERPRINT_Q, "_featureSpace_pairwise_tests.csv")))

plot_feature <- make_dist2gr_plot(
  dist_df = dist_feature,
  focus_pairs = FOCUS_PAIRS,
  ylab = "Distance to GR centroid (AB signature feature space)"
)

ggsave(file.path(OUTDIR, "05_Dist2GR", paste0("Dist2GR_ABsig_q", FINGERPRINT_Q, "_featureSpace.png")),
       plot_feature$plot, width = 7.6, height = 5.2, dpi = 300)
ggsave(file.path(OUTDIR, "05_Dist2GR", paste0("Dist2GR_ABsig_q", FINGERPRINT_Q, "_featureSpace.pdf")),
       plot_feature$plot, width = 7.6, height = 5.2)

# ----------------------------
# 7) Session info
# ----------------------------
sink(file.path(OUTDIR, "sessionInfo.txt"))
print(sessionInfo())
sink()

msg("Done.")
msg("Output directory: ", normalizePath(OUTDIR, winslash = "/"))