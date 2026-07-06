#!/usr/bin/env Rscript
# =========================================================
# 3-group targeted metabolomics | stats + volcano + heatmap
# Sankey-aligned statistical framework
# ---------------------------------------------------------
# Input:
#   - Q2.xlsx (sheet = "sample")
# Groups:
#   A = GR
#   B = MR
#   C = MR+Phage cocktail
# Statistical framework (aligned to current Sankey/Venn line):
#   * 0 treated as not detected and imputed to half of the minimum
#     positive value PER metabolite across all retained 3-group samples
#   * log2 transform
#   * overall one-way ANOVA across 3 groups -> p_overall; BH -> fdr_overall
#   * pairwise pooled-LSD using ANOVA pooled MSE -> p_* ; BH within contrast -> fdr_*
#   * log2FC = difference of group log2 means
#   * gm_* = 2^(mean_*)
#   * rescue_score_C_toward_A = closeness of C toward A relative to A vs B distance
# Outputs:
#   results_3group_metabolomics_stats_volcano_heatmap_sankey_aligned/
#     01_stats/
#     02_volcano/
#     03_heatmap/
# =========================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(pheatmap)
})

HAS_GGREPEL <- requireNamespace("ggrepel", quietly = TRUE)
if (HAS_GGREPEL) suppressPackageStartupMessages(library(ggrepel))
HAS_OPENXLSX <- requireNamespace("openxlsx", quietly = TRUE)

# ----------------------------
# User parameters
# ----------------------------
MET_QUANT_XLSX <- "Q2.xlsx"
SHEET_NAME <- "sample"

OUTDIR <- "results_3group_metabolomics_stats_volcano_heatmap_sankey_aligned"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "01_stats"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "02_volcano"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "03_heatmap"), showWarnings = FALSE, recursive = TRUE)

# Signature / plotting threshold (same as your previous metabolite signature line)
AB_SIGNATURE_Q_CUT <- 0.20
VOLCANO_FC_LINE <- 1.0  # reference line only

# Heatmap options
HEATMAP_CLUSTER_COLS <- FALSE
HEATMAP_SHOW_SAMPLE_NAMES <- FALSE
HEATMAP_ROW_ORDER_BY_LOG2FC <- TRUE
HEATMAP_USE_CHINESE_LABEL <- FALSE  # safer for cross-platform plotting

GROUP_LEVELS <- c("GR", "MR", "MR+Phage cocktail")
GROUP_COLORS <- c(
  "GR" = "#7CAE00",
  "MR" = "#F8766D",
  "MR+Phage cocktail" = "#00BFC4"
)

# ----------------------------
# Helpers
# ----------------------------
msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n")

lsd_p <- function(mean1, mean2, n1, n2, mse, df_res) {
  if (any(is.na(c(mean1, mean2, n1, n2, mse, df_res))) ||
      n1 < 2 || n2 < 2 || mse <= 0 || df_res <= 0) {
    return(NA_real_)
  }
  tval <- (mean1 - mean2) / sqrt(mse * (1 / n1 + 1 / n2))
  2 * stats::pt(abs(tval), df = df_res, lower.tail = FALSE)
}

safe_write_csv <- function(x, file) {
  readr::write_csv(x, file, na = "")
}

rescue_score <- function(mean_A, mean_B, mean_C) {
  denom <- abs(mean_A - mean_B)
  if (is.na(denom) || denom == 0) return(NA_real_)
  (abs(mean_A - mean_B) - abs(mean_A - mean_C)) / denom
}

impute_half_minpos <- function(x) {
  x <- as.numeric(x)
  min_pos <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
  if (!is.finite(min_pos)) return(rep(NA_real_, length(x)))
  x[is.na(x) | x == 0] <- min_pos / 2
  x
}

row_z <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sdv <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) return(rep(0, length(x)))
  (x - mu) / sdv
}

pick_label <- function(metabolite, zh_name, use_chinese = FALSE) {
  if (use_chinese && !is.na(zh_name) && zh_name != "") return(as.character(zh_name))
  as.character(metabolite)
}

# ----------------------------
# Read data
# ----------------------------
if (!file.exists(MET_QUANT_XLSX)) stop("File not found: ", MET_QUANT_XLSX)

msg("Reading metabolomics raw table...")
df_raw <- readxl::read_excel(MET_QUANT_XLSX, sheet = SHEET_NAME)
df_raw <- as.data.frame(df_raw)

sample_cols <- grep("^[ABC][1-9][0-9]*$", names(df_raw), value = TRUE)
if (length(sample_cols) == 0) stop("Cannot find sample columns like A1..C5 in ", MET_QUANT_XLSX)

# stable A->B->C order
sample_cols <- c(sample_cols[grepl("^A", sample_cols)],
                 sample_cols[grepl("^B", sample_cols)],
                 sample_cols[grepl("^C", sample_cols)])

need_cols <- c("Compound Name (nmol/g)", sample_cols)
miss <- setdiff(need_cols, names(df_raw))
if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))

meta_cols <- intersect(
  c("中文名称", "CAS", "Formula", "Pubchem.ID", "HMDB", "KEGG", "Class", "分类", "Average Molecular Weight", "average Molecular'"),
  names(df_raw)
)

sample_meta <- tibble::tibble(
  sample = sample_cols,
  group = dplyr::case_when(
    grepl("^A", sample) ~ "GR",
    grepl("^B", sample) ~ "MR",
    grepl("^C", sample) ~ "MR+Phage cocktail",
    TRUE ~ NA_character_
  )
) %>%
  filter(!is.na(group)) %>%
  mutate(group = factor(group, levels = GROUP_LEVELS))

safe_write_csv(sample_meta, file.path(OUTDIR, "01_stats", "sample_metadata_used.csv"))

# ----------------------------
# Long table + imputation
# ----------------------------
met_long <- df_raw %>%
  select(all_of(c("Compound Name (nmol/g)", meta_cols, sample_cols))) %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "sample", values_to = "value_raw") %>%
  mutate(
    metabolite = as.character(`Compound Name (nmol/g)`),
    value_raw = suppressWarnings(as.numeric(value_raw))
  ) %>%
  left_join(sample_meta, by = "sample") %>%
  filter(!is.na(group), !is.na(metabolite), metabolite != "")

half_min_tbl <- met_long %>%
  group_by(metabolite) %>%
  summarise(
    min_pos = suppressWarnings(min(value_raw[value_raw > 0], na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    min_pos = ifelse(is.infinite(min_pos), NA_real_, min_pos),
    half_min = min_pos / 2
  )

met_long2 <- met_long %>%
  left_join(half_min_tbl, by = "metabolite") %>%
  mutate(
    value_imp = case_when(
      is.na(value_raw) ~ NA_real_,
      value_raw > 0 ~ value_raw,
      TRUE ~ half_min
    ),
    value_log2 = ifelse(!is.na(value_imp) & value_imp > 0, log2(value_imp), NA_real_)
  )

# ----------------------------
# Statistics table (3-group, Sankey-aligned)
# ----------------------------
msg("Computing full statistics table...")

stats_tbl <- met_long2 %>%
  group_by(metabolite) %>%
  group_modify(function(d, key) {
    gsum <- d %>%
      group_by(group) %>%
      summarise(
        n = sum(!is.na(value_log2)),
        mean = mean(value_log2, na.rm = TRUE),
        sd = sd(value_log2, na.rm = TRUE),
        .groups = "drop"
      )

    get_g <- function(g, field) {
      x <- gsum %>% filter(group == g) %>% pull({{ field }})
      if (length(x) == 0) NA_real_ else x[1]
    }

    mean_A <- get_g("GR", mean); sd_A <- get_g("GR", sd); nA <- get_g("GR", n)
    mean_B <- get_g("MR", mean); sd_B <- get_g("MR", sd); nB <- get_g("MR", n)
    mean_C <- get_g("MR+Phage cocktail", mean); sd_C <- get_g("MR+Phage cocktail", sd); nC <- get_g("MR+Phage cocktail", n)

    p_overall <- NA_real_
    mse <- NA_real_
    df_res <- NA_real_
    dd <- d %>% filter(!is.na(value_log2))
    if (nrow(dd) >= 6 && length(unique(dd$group)) >= 2) {
      fit <- try(stats::aov(value_log2 ~ group, data = dd), silent = TRUE)
      if (!inherits(fit, "try-error")) {
        s <- summary(fit)[[1]]
        p_overall <- s[["Pr(>F)"]][1]
        mse <- s[["Mean Sq"]][nrow(s)]
        df_res <- s[["Df"]][nrow(s)]
      }
    }

    p_A_vs_B <- lsd_p(mean_A, mean_B, nA, nB, mse, df_res)
    p_C_vs_B <- lsd_p(mean_C, mean_B, nC, nB, mse, df_res)
    p_C_vs_A <- lsd_p(mean_C, mean_A, nC, nA, mse, df_res)

    log2FC_A_vs_B <- mean_A - mean_B
    log2FC_C_vs_B <- mean_C - mean_B
    log2FC_C_vs_A <- mean_C - mean_A

    gm_A <- 2 ^ mean_A
    gm_B <- 2 ^ mean_B
    gm_C <- 2 ^ mean_C

    rescue_C_toward_A <- rescue_score(mean_A, mean_B, mean_C)

    tibble::tibble(
      p_overall = p_overall,
      mean_A = mean_A, sd_A = sd_A,
      mean_B = mean_B, sd_B = sd_B,
      mean_C = mean_C, sd_C = sd_C,
      log2FC_A_vs_B = log2FC_A_vs_B, p_A_vs_B = p_A_vs_B,
      log2FC_C_vs_B = log2FC_C_vs_B, p_C_vs_B = p_C_vs_B,
      log2FC_C_vs_A = log2FC_C_vs_A, p_C_vs_A = p_C_vs_A,
      gm_A = gm_A, gm_B = gm_B, gm_C = gm_C,
      rescue_score_C_toward_A = rescue_C_toward_A
    )
  }) %>%
  ungroup()

meta_one <- df_raw %>%
  transmute(
    metabolite = as.character(`Compound Name (nmol/g)`),
    中文名称 = if ("中文名称" %in% names(df_raw)) as.character(`中文名称`) else NA_character_,
    CAS = if ("CAS" %in% names(df_raw)) as.character(CAS) else NA_character_,
    formula = if ("Formula" %in% names(df_raw)) as.character(Formula) else NA_character_,
    pubchemid = if ("Pubchem.ID" %in% names(df_raw)) as.character(`Pubchem.ID`) else NA_character_,
    hmdb = if ("HMDB" %in% names(df_raw)) as.character(HMDB) else NA_character_,
    kegg = if ("KEGG" %in% names(df_raw)) as.character(KEGG) else NA_character_,
    class = if ("Class" %in% names(df_raw)) as.character(Class) else NA_character_,
    分类 = if ("分类" %in% names(df_raw)) as.character(`分类`) else NA_character_,
    average_molecular_weight = if ("Average Molecular Weight" %in% names(df_raw)) suppressWarnings(as.numeric(`Average Molecular Weight`)) else if ("average Molecular'" %in% names(df_raw)) suppressWarnings(as.numeric(`average Molecular'`)) else NA_real_
  ) %>%
  distinct()

out <- meta_one %>%
  left_join(stats_tbl, by = "metabolite") %>%
  mutate(
    fdr_overall = p.adjust(p_overall, method = "BH"),
    fdr_A_vs_B = p.adjust(p_A_vs_B, method = "BH"),
    fdr_C_vs_B = p.adjust(p_C_vs_B, method = "BH"),
    fdr_C_vs_A = p.adjust(p_C_vs_A, method = "BH")
  ) %>%
  select(
    metabolite,
    p_overall,
    mean_A, sd_A, mean_B, sd_B, mean_C, sd_C,
    log2FC_A_vs_B, p_A_vs_B,
    log2FC_C_vs_B, p_C_vs_B,
    log2FC_C_vs_A, p_C_vs_A,
    fdr_A_vs_B, fdr_C_vs_B, fdr_C_vs_A, fdr_overall,
    CAS, formula, pubchemid, hmdb, kegg, class,
    average_molecular_weight,
    gm_A, gm_B, gm_C,
    rescue_score_C_toward_A,
    中文名称, 分类
  )

safe_write_csv(out, file.path(OUTDIR, "01_stats", "metabolomics_stats_all_results_3groups_sankey_aligned.csv"))
if (HAS_OPENXLSX) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "all_results")
  openxlsx::writeData(wb, "all_results", out)
  openxlsx::saveWorkbook(wb, file.path(OUTDIR, "01_stats", "metabolomics_stats_all_results_3groups_sankey_aligned.xlsx"), overwrite = TRUE)
}

# ----------------------------
# Volcano plot (GR vs MR)
# ----------------------------
msg("Drawing volcano plot...")
vol_dat <- out %>%
  transmute(
    metabolite = as.character(metabolite),
    label = vapply(seq_len(n()), function(i) pick_label(metabolite[i], 中文名称[i], use_chinese = FALSE), character(1)),
    log2FC = as.numeric(log2FC_A_vs_B),
    pval = as.numeric(p_A_vs_B),
    qval = as.numeric(fdr_A_vs_B)
  ) %>%
  filter(is.finite(log2FC), is.finite(pval), is.finite(qval), pval > 0, qval > 0) %>%
  mutate(
    is_sig = qval < AB_SIGNATURE_Q_CUT,
    class = case_when(
      is_sig & log2FC > 0 ~ "Up (q<0.2)",
      is_sig & log2FC < 0 ~ "Down (q<0.2)",
      TRUE ~ "Not significant"
    ),
    y = -log10(qval),
    label_sig = ifelse(is_sig, label, NA_character_)
  )

vol_dat$class <- factor(vol_dat$class, levels = c("Down (q<0.2)", "Not significant", "Up (q<0.2)"))

p_vol <- ggplot(vol_dat, aes(x = log2FC, y = y, color = class)) +
  geom_point(aes(alpha = class), size = 2.1) +
  scale_alpha_manual(values = c(
    "Down (q<0.2)" = 0.95,
    "Not significant" = 0.75,
    "Up (q<0.2)" = 0.95
  ), guide = "none") +
  scale_color_manual(values = c(
    "Up (q<0.2)" = "#E41A1C",
    "Down (q<0.2)" = "#377EB8",
    "Not significant" = "grey50"
  )) +
  geom_hline(yintercept = -log10(AB_SIGNATURE_Q_CUT), linetype = 2, linewidth = 0.7, color = "grey50") +
  geom_vline(xintercept = c(-VOLCANO_FC_LINE, VOLCANO_FC_LINE), linetype = 2, linewidth = 0.7, color = "grey50") +
  labs(
    x = "log2FC (GR vs MR)",
    y = expression(-log[10](q)),
    color = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = c(0.98, 0.85),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA)
  ) +
  guides(
    color = guide_legend(
      override.aes = list(alpha = 1, shape = 16, size = 4, stroke = 0)
    )
  )

sig_dat <- subset(vol_dat, !is.na(label_sig))
if (nrow(sig_dat) > 0) {
  if (HAS_GGREPEL) {
    p_vol <- p_vol + ggrepel::geom_text_repel(
      data = sig_dat,
      aes(label = label_sig),
      size = 4.0,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.25,
      segment.color = "grey40",
      min.segment.length = 0
    )
  } else {
    p_vol <- p_vol + geom_text(data = sig_dat, aes(label = label_sig), size = 4.0, vjust = -0.6)
  }
}

ggsave(file.path(OUTDIR, "02_volcano", "Volcano_GR_vs_MR_white_yQ_3groups.png"), p_vol, width = 10, height = 7, dpi = 300, bg = "white")
ggsave(file.path(OUTDIR, "02_volcano", "Volcano_GR_vs_MR_white_yQ_3groups.pdf"), p_vol, width = 10, height = 7, bg = "white")

# ----------------------------
# Heatmap (AB signature, 3 groups)
# ----------------------------
msg("Drawing AB-signature heatmap...")
sig_tbl <- out %>%
  transmute(
    metabolite = as.character(metabolite),
    q = as.numeric(fdr_A_vs_B),
    log2FC = as.numeric(log2FC_A_vs_B),
    中文名称 = as.character(中文名称)
  ) %>%
  filter(is.finite(q), is.finite(log2FC), q < AB_SIGNATURE_Q_CUT) %>%
  mutate(
    row_label = mapply(pick_label, metabolite, 中文名称,
                       MoreArgs = list(use_chinese = HEATMAP_USE_CHINESE_LABEL),
                       SIMPLIFY = TRUE, USE.NAMES = FALSE)
  ) %>%
  arrange(log2FC)

if (nrow(sig_tbl) == 0) {
  warning("No AB signature metabolites found at q<0.2. Heatmap skipped.")
} else {
  dat0 <- df_raw %>%
    select(all_of(c("Compound Name (nmol/g)", meta_cols, sample_cols))) %>%
    mutate(across(all_of(sample_cols), ~ suppressWarnings(as.numeric(.x))))

  met_names <- as.character(dat0[["Compound Name (nmol/g)"]])
  mat <- as.matrix(dat0[, sample_cols, drop = FALSE])
  rownames(mat) <- make.unique(met_names)

  keep <- sig_tbl$metabolite[sig_tbl$metabolite %in% rownames(mat)]
  if (length(keep) == 0) {
    warning("Signature metabolites not found in quantitative table. Heatmap skipped.")
  } else {
    mat <- mat[keep, , drop = FALSE]

    mat_imp <- t(apply(mat, 1, impute_half_minpos))
    colnames(mat_imp) <- colnames(mat)
    rownames(mat_imp) <- rownames(mat)

    mat_log2 <- log2(mat_imp)
    mat_z <- t(apply(mat_log2, 1, row_z))
    colnames(mat_z) <- colnames(mat_log2)
    rownames(mat_z) <- rownames(mat_log2)

    sample_order <- c(sample_cols[grepl("^A", sample_cols)],
                      sample_cols[grepl("^B", sample_cols)],
                      sample_cols[grepl("^C", sample_cols)])
    sample_order <- sample_order[sample_order %in% colnames(mat_z)]
    mat_z <- mat_z[, sample_order, drop = FALSE]

    group_raw <- substr(colnames(mat_z), 1, 1)
    group_lab <- dplyr::recode(
      group_raw,
      "A" = "GR",
      "B" = "MR",
      "C" = "MR+Phage cocktail"
    )

    anno_col <- data.frame(Group = factor(group_lab, levels = GROUP_LEVELS))
    rownames(anno_col) <- colnames(mat_z)

    gaps_col <- c(sum(grepl("^A", colnames(mat_z))), sum(grepl("^[AB]", colnames(mat_z))))
    gaps_col <- gaps_col[gaps_col < ncol(mat_z)]

    cluster_rows_flag <- !HEATMAP_ROW_ORDER_BY_LOG2FC

    if (HEATMAP_ROW_ORDER_BY_LOG2FC) {
      ord_tbl <- sig_tbl %>% filter(metabolite %in% rownames(mat_z)) %>% arrange(log2FC)
      mat_z <- mat_z[ord_tbl$metabolite, , drop = FALSE]
      rownames(mat_z) <- ord_tbl$row_label
    }

    pheatmap(
      mat_z,
      annotation_col = anno_col,
      annotation_colors = list(Group = GROUP_COLORS),
      cluster_rows = cluster_rows_flag,
      cluster_cols = HEATMAP_CLUSTER_COLS,
      gaps_col = if (!HEATMAP_CLUSTER_COLS) gaps_col else NULL,
      show_colnames = HEATMAP_SHOW_SAMPLE_NAMES,
      show_rownames = TRUE,
      fontsize_row = 10,
      fontsize_col = 9,
      border_color = NA,
      filename = file.path(OUTDIR, "03_heatmap", "Heatmap_ABsignature_q0.2_3groups.png"),
      width = 10,
      height = 8
    )

    pheatmap(
      mat_z,
      annotation_col = anno_col,
      annotation_colors = list(Group = GROUP_COLORS),
      cluster_rows = cluster_rows_flag,
      cluster_cols = HEATMAP_CLUSTER_COLS,
      gaps_col = if (!HEATMAP_CLUSTER_COLS) gaps_col else NULL,
      show_colnames = HEATMAP_SHOW_SAMPLE_NAMES,
      show_rownames = TRUE,
      fontsize_row = 10,
      fontsize_col = 9,
      border_color = NA,
      filename = file.path(OUTDIR, "03_heatmap", "Heatmap_ABsignature_q0.2_3groups.pdf"),
      width = 10,
      height = 8
    )

    heatmap_per_sample <- as.data.frame(mat_log2[keep, sample_order, drop = FALSE]) %>%
      tibble::rownames_to_column("metabolite") %>%
      left_join(sig_tbl %>% select(metabolite, row_label, q, log2FC), by = "metabolite")
    safe_write_csv(heatmap_per_sample, file.path(OUTDIR, "03_heatmap", "Heatmap_ABsignature_q0.2_3groups_input_log2.csv"))
  }
}

msg("Done.")
msg("Output directory: ", normalizePath(OUTDIR, winslash = "/"))
