# ============================================================
# 02_PCA_export_tables_metabolomics_FULL6.R  (NO openxlsx)
# PCA exports (All metabolites + Fingerprint metabolites) with:
# - sample scores (PC1..PC5)
# - variance explained
# - loadings
# - pairwise full 6 comparisons on PC scores (Wilcoxon): stat + p + BH-q
# - distance-to-GR-centroid (PC1-2, PC1-5) per sample + pairwise full 6
#
# Inputs (in working directory):
#   1) 定量结果2.xlsx   (sheet: 样本浓度)
#   2) neurotransmitter_metabolomics_stats_q0_2.xlsx (sheet: all_results)
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(readr)
})

# =========================
# 0) PARAMETERS
# =========================
MET_QUANT_XLSX  <- "Q2.xlsx"
MET_QUANT_SHEET <- "sample"

STATS_XLSX  <- "neurotransmitter_metabolomics_stats_q0_2.xlsx"
STATS_SHEET <- "all_results"

FINGERPRINT_Q <- 0.2
PC_N <- 5

OUTDIR <- "Metabolomics_PCA_exports_FULL6"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Group mapping from sample IDs A1..D5
GROUP_LEVELS <- c("GR","MR","MR+Phage cocktail","MR+AFMT")
map_group <- function(sample_id){
  g <- substr(sample_id, 1, 1)
  dplyr::recode(g,
         "A"="GR",
         "B"="MR",
         "C"="MR+Phage cocktail",
         "D"="MR+AFMT",
         .default = NA_character_)
}

# =========================
# 1) HELPERS
# =========================
pick_col_regex <- function(nms, patterns) {
  for (pat in patterns) {
    hit <- grep(pat, nms, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

impute_row_halfmin <- function(v){
  v <- suppressWarnings(as.numeric(v))
  if (all(is.na(v))) return(rep(NA_real_, length(v)))
  min_pos <- suppressWarnings(min(v[v > 0], na.rm = TRUE))
  if (is.infinite(min_pos) || is.na(min_pos)) return(rep(NA_real_, length(v)))
  v[v == 0] <- min_pos / 2
  v
}

run_pca_from_quant <- function(raw_df, sample_cols, feature_keep = NULL, feature_name_col = "metabolite", pc_n = 5) {

  df <- raw_df %>%
    mutate(metabolite = as.character(.data[[feature_name_col]])) %>%
    select(metabolite, all_of(sample_cols))

  if (!is.null(feature_keep)) {
    df <- df %>% filter(metabolite %in% feature_keep)
  }

  # impute zeros per row, log2
  df2 <- df %>%
    rowwise() %>%
    mutate(across(all_of(sample_cols), ~ impute_row_halfmin(.x))) %>%
    ungroup() %>%
    mutate(across(all_of(sample_cols), ~ log2(.x)))

  # matrix rows=metabolites, cols=samples
  mat <- df2 %>% as.data.frame()
  rownames(mat) <- make.unique(mat$metabolite)
  mat$metabolite <- NULL
  X <- as.matrix(mat)

  # keep finite & variable
  keep <- apply(X, 1, function(x) all(is.finite(x)) && sd(x) > 0)
  X <- X[keep, , drop = FALSE]

  # row z-score
  Xz <- t(scale(t(X)))
  Xz[!is.finite(Xz)] <- 0

  # PCA on samples
  pca <- prcomp(t(Xz), center = FALSE, scale. = FALSE)

  k <- min(pc_n, ncol(pca$x))
  scores <- as.data.frame(pca$x[, 1:k, drop = FALSE]) %>%
    rownames_to_column("sample") %>%
    mutate(group = factor(map_group(sample), levels = GROUP_LEVELS))

  eig <- pca$sdev^2
  var_pct <- eig / sum(eig) * 100
  variance <- tibble(
    PC = paste0("PC", seq_along(eig)),
    sdev = as.numeric(pca$sdev),
    eigenvalue = as.numeric(eig),
    variance_pct = as.numeric(var_pct)
  )

  rot <- as.data.frame(pca$rotation[, 1:k, drop = FALSE]) %>%
    rownames_to_column("metabolite")

  list(scores = scores, variance = variance, loadings = rot, Xz = Xz)
}

pairwise_wilcox_full6 <- function(df, value_col, group_col = "group") {
  df <- df %>% filter(!is.na(.data[[value_col]])) %>%
    mutate(group = droplevels(as.factor(.data[[group_col]])))

  grps <- levels(df$group)
  if (length(grps) != 4) stop("Expected 4 groups for full 6 comparisons; got: ", paste(grps, collapse = ", "))

  pairs <- t(combn(grps, 2))

  out <- bind_rows(apply(pairs, 1, function(x) {
    g1 <- x[1]; g2 <- x[2]
    v1 <- df %>% filter(group == g1) %>% pull(.data[[value_col]])
    v2 <- df %>% filter(group == g2) %>% pull(.data[[value_col]])
    wt <- suppressWarnings(wilcox.test(v1, v2, exact = FALSE))

    tibble(
      group1 = g1, group2 = g2,
      n1 = length(v1), n2 = length(v2),
      median1 = median(v1, na.rm = TRUE),
      median2 = median(v2, na.rm = TRUE),
      delta_median = median1 - median2,
      stat_W = unname(wt$statistic),
      p_value = wt$p.value
    )
  })) %>%
    mutate(q_value = p.adjust(p_value, method = "BH"))

  out
}

dist_to_GR_centroid <- function(scores_df, pcs = c("PC1","PC2")) {
  stopifnot(all(pcs %in% names(scores_df)))
  gr_cent <- scores_df %>% filter(group == "GR") %>%
    summarise(across(all_of(pcs), mean, na.rm = TRUE))

  center <- as.numeric(gr_cent[1, pcs, drop = TRUE])
  mat <- as.matrix(scores_df[, pcs, drop = FALSE])
  d <- sqrt(rowSums((mat - matrix(center, nrow = nrow(mat), ncol = length(center), byrow = TRUE))^2))

  tibble(
    sample = scores_df$sample,
    group  = scores_df$group,
    pcs_used = paste(pcs, collapse = "+"),
    dist_to_GR = as.numeric(d)
  )
}

# =========================
# 2) READ INPUTS
# =========================
stopifnot(file.exists(MET_QUANT_XLSX))
raw <- readxl::read_excel(MET_QUANT_XLSX, sheet = MET_QUANT_SHEET)

sample_cols <- names(raw)[str_detect(names(raw), "^[ABCD]\\d+$")]
if (length(sample_cols) == 0) stop("Cannot detect sample columns like A1..D5 in quant table.")

met_col_quant <- pick_col_regex(
  names(raw),
  patterns = c("^Compound Name \\(nmol/g\\)$", "^Compound Name", "^metabolite$", "^Metabolite$", "^Name$", "^代谢物$")
)
if (is.na(met_col_quant)) met_col_quant <- names(raw)[1]

raw2 <- raw %>%
  mutate(metabolite = as.character(.data[[met_col_quant]])) %>%
  select(metabolite, all_of(sample_cols))

stopifnot(file.exists(STATS_XLSX))
stats <- readxl::read_excel(STATS_XLSX, sheet = STATS_SHEET)

q_col_ab <- pick_col_regex(
  names(stats),
  patterns = c("^fdr_a_vs_b$", "^fdr_a_v__b$", "fdr.*a.*b")
)
if (is.na(q_col_ab)) stop("Cannot find GR vs MR FDR column in stats. Please check column names.")

met_col_stats <- pick_col_regex(names(stats), patterns = c("^metabolite$","^Metabolite$","^name$"))
if (is.na(met_col_stats)) stop("Cannot find metabolite name column in stats.")

fingerprint_mets <- stats %>%
  transmute(metabolite = as.character(.data[[met_col_stats]]),
            q_ab = suppressWarnings(as.numeric(.data[[q_col_ab]]))) %>%
  filter(!is.na(metabolite) & metabolite != "", !is.na(q_ab), q_ab < FINGERPRINT_Q) %>%
  distinct(metabolite) %>%
  pull(metabolite)

# Save fingerprint list
write_csv(tibble(metabolite = fingerprint_mets),
          file.path(OUTDIR, paste0("Fingerprint_list_q", FINGERPRINT_Q, ".csv")))

# =========================
# 3) RUN PCA
# =========================
res_all <- run_pca_from_quant(raw2, sample_cols, feature_keep = NULL, feature_name_col = "metabolite", pc_n = PC_N)
res_fp  <- run_pca_from_quant(raw2, sample_cols, feature_keep = fingerprint_mets, feature_name_col = "metabolite", pc_n = PC_N)

# =========================
# 4) EXPORT
# =========================
export_one_pca <- function(res, tag) {
  write_csv(res$scores,   file.path(OUTDIR, paste0("PCA_", tag, "_scores_PC1toPC", PC_N, ".csv")))
  write_csv(res$variance, file.path(OUTDIR, paste0("PCA_", tag, "_variance.csv")))
  write_csv(res$loadings, file.path(OUTDIR, paste0("PCA_", tag, "_loadings_PC1toPC", PC_N, ".csv")))

  pcs <- paste0("PC", 1:min(PC_N, ncol(res$scores) - 2))

  cent <- res$scores %>%
    group_by(group) %>%
    summarise(across(all_of(pcs), mean, na.rm = TRUE), .groups = "drop")
  write_csv(cent, file.path(OUTDIR, paste0("PCA_", tag, "_group_centroids.csv")))

  # pairwise tests per PC (full 6; BH within each PC)
  pw_pc <- lapply(pcs, function(pc) {
    pairwise_wilcox_full6(res$scores, value_col = pc, group_col = "group") %>%
      mutate(dataset = tag, PC = pc) %>%
      select(dataset, PC, everything())
  }) %>% bind_rows()
  write_csv(pw_pc, file.path(OUTDIR, paste0("PCA_", tag, "_pairwise_Wilcoxon_full6_PC_scores.csv")))

  # distances to GR centroid
  pcs12 <- intersect(c("PC1","PC2"), pcs)
  d12 <- dist_to_GR_centroid(res$scores, pcs = pcs12)
  d15 <- dist_to_GR_centroid(res$scores, pcs = pcs)

  write_csv(d12, file.path(OUTDIR, paste0("PCA_", tag, "_dist_to_GR_PC1_2_perSample.csv")))
  write_csv(d15, file.path(OUTDIR, paste0("PCA_", tag, "_dist_to_GR_PC1_5_perSample.csv")))

  pw_d12 <- pairwise_wilcox_full6(d12, value_col = "dist_to_GR", group_col = "group") %>%
    mutate(dataset = tag, pcs_used = "PC1+PC2") %>% select(dataset, pcs_used, everything())
  pw_d15 <- pairwise_wilcox_full6(d15, value_col = "dist_to_GR", group_col = "group") %>%
    mutate(dataset = tag, pcs_used = paste(pcs, collapse = "+")) %>% select(dataset, pcs_used, everything())

  write_csv(pw_d12, file.path(OUTDIR, paste0("PCA_", tag, "_pairwise_Wilcoxon_full6_dist_to_GR_PC1_2.csv")))
  write_csv(pw_d15, file.path(OUTDIR, paste0("PCA_", tag, "_pairwise_Wilcoxon_full6_dist_to_GR_PC1_5.csv")))

  invisible(TRUE)
}

export_one_pca(res_all, "AllMet")
export_one_pca(res_fp,  paste0("Fingerprint_q", FINGERPRINT_Q))

sink(file.path(OUTDIR, "sessionInfo.txt"))
print(sessionInfo())
sink()

cat("[DONE] PCA export tables saved to:", OUTDIR, "\n")
cat("Fingerprint size (q<", FINGERPRINT_Q, "): ", length(fingerprint_mets), "\n", sep = "")
