#!/usr/bin/env Rscript
# =========================================================
# Targeted metabolomics (AQ700) | Full statistics table (platform-aligned)
# - Input: "定量结果2.xlsx" (sheet: 样本浓度)
# - Output: neurotransmitter_metabolomics_stats_reproduced.xlsx
# - Groups: A=GR, B=MR, C=MR+Phage cocktail, D=MR+AFMT
# - Statistics (platform口径):
#   * zero (0) treated as "not detected" and imputed to half of the minimum
#     positive value PER metabolite across all samples (LOD/2 proxy) before log2
#   * log2 transform
#   * overall one-way ANOVA across 4 groups -> p_overall; BH -> fdr_overall
#   * pairwise Fisher's LSD using ANOVA pooled MSE (common variance) -> p_*
#     BH across metabolites within each contrast -> fdr_*
#   * log2FC = difference of group log2 means (equivalent to log2 ratio of geometric means)
#   * gm_* = 2^(mean_*)
#   * rescue_% computed on geometric means relative to GR vs MR baseline
# =========================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(openxlsx)
})

# ----------------------------
# User parameters
# ----------------------------
MET_QUANT_XLSX <- "定量结果2.xlsx"
SHEET_NAME <- "样本浓度"

OUT_XLSX <- "neurotransmitter_metabolomics_stats_reproduced.xlsx"
OUT_SHEET <- "all_results"

# Sample columns expected (A1-5, B1-5, C1-5, D1-5)
sample_cols <- c(paste0("A", 1:5), paste0("B", 1:5), paste0("C", 1:5), paste0("D", 1:5))

# ----------------------------
# Helpers
# ----------------------------
stop_if_missing <- function(x, what="required columns") {
  miss <- setdiff(x, names(df_raw))
  if (length(miss) > 0) stop("Missing ", what, ": ", paste(miss, collapse=", "))
}

# Fisher LSD p-value for mean difference using pooled ANOVA MSE
lsd_p <- function(mean1, mean2, n1, n2, mse, df_res) {
  if (any(is.na(c(mean1, mean2, n1, n2, mse, df_res))) || n1 < 2 || n2 < 2 || mse <= 0 || df_res <= 0) return(NA_real_)
  tval <- (mean1 - mean2) / sqrt(mse * (1/n1 + 1/n2))
  2 * stats::pt(abs(tval), df = df_res, lower.tail = FALSE)
}

# Rescue % on geometric means (gm); direction determined by GR vs MR
rescue_pct <- function(gm_gr, gm_mr, gm_treat) {
  if (any(is.na(c(gm_gr, gm_mr, gm_treat))) || (gm_gr == gm_mr)) return(NA_real_)
  if (gm_mr > gm_gr) { # disease-up, rescue means decrease back to GR
    (gm_mr - gm_treat) / (gm_mr - gm_gr) * 100
  } else {            # disease-down, rescue means increase back to GR
    (gm_treat - gm_mr) / (gm_gr - gm_mr) * 100
  }
}

# Safe numeric conversion
as_num <- function(x) suppressWarnings(as.numeric(x))

# ----------------------------
# Read data
# ----------------------------
if (!file.exists(MET_QUANT_XLSX)) stop("File not found: ", MET_QUANT_XLSX)

sheets <- readxl::excel_sheets(MET_QUANT_XLSX)
if (!(SHEET_NAME %in% sheets)) {
  message("[WARN] Sheet '", SHEET_NAME, "' not found. Using first sheet: ", sheets[1])
  SHEET_NAME <- sheets[1]
}

df_raw <- readxl::read_excel(MET_QUANT_XLSX, sheet = SHEET_NAME)
df_raw <- as.data.frame(df_raw)

# Basic sanity
need_cols <- c("Compound Name (nmol/g)", sample_cols)
miss <- setdiff(need_cols, names(df_raw))
if (length(miss) > 0) stop("Missing required columns in quant file: ", paste(miss, collapse=", "))

# Keep meta columns if present (we do not force them)
meta_cols <- intersect(c("中文名称","CAS","Formula","Pubchem.ID","HMDB","KEGG","Class","分类","average Molecular'"), names(df_raw))

# Long format
df_long <- df_raw %>%
  dplyr::select(all_of(c("Compound Name (nmol/g)", meta_cols, sample_cols))) %>%
  tidyr::pivot_longer(cols = all_of(sample_cols), names_to = "sample", values_to = "value_raw") %>%
  dplyr::mutate(
    metabolite = as.character(`Compound Name (nmol/g)`),
    value_raw = as_num(value_raw),
    group_code = substr(sample, 1, 1),
    group = dplyr::recode(group_code,
                          "A" = "GR",
                          "B" = "MR",
                          "C" = "MR_Phage",
                          "D" = "MR_AFMT",
                          .default = NA_character_)
  ) %>%
  dplyr::filter(!is.na(group) & !is.na(metabolite) & metabolite != "")

# ----------------------------
# Per-metabolite imputation + stats
# ----------------------------
# Compute per-metabolite half-min positive for imputation (LOD/2 proxy)
half_min_pos <- df_long %>%
  group_by(metabolite) %>%
  summarize(min_pos = suppressWarnings(min(value_raw[value_raw > 0], na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(min_pos = ifelse(is.infinite(min_pos), NA_real_, min_pos),
         half_min = min_pos / 2)

df_long2 <- df_long %>%
  left_join(half_min_pos, by = "metabolite") %>%
  mutate(
    value_imp = case_when(
      is.na(value_raw) ~ NA_real_,
      value_raw > 0 ~ value_raw,
      TRUE ~ half_min
    ),
    value_log2 = ifelse(!is.na(value_imp) & value_imp > 0, log2(value_imp), NA_real_)
  )

# Core stats per metabolite
stats_tbl <- df_long2 %>%
  group_by(metabolite) %>%
  group_modify(function(d, key) {
    # group summaries on log2 scale
    gsum <- d %>%
      group_by(group) %>%
      summarize(
        n = sum(!is.na(value_log2)),
        mean = mean(value_log2, na.rm = TRUE),
        sd = sd(value_log2, na.rm = TRUE),
        .groups = "drop"
      )

    # Ensure all groups present
    get_g <- function(g, field) {
      x <- gsum %>% filter(group == g) %>% pull({{field}})
      if (length(x) == 0) NA_real_ else x[1]
    }

    mean_A <- get_g("GR", mean); sd_A <- get_g("GR", sd); nA <- get_g("GR", n)
    mean_B <- get_g("MR", mean); sd_B <- get_g("MR", sd); nB <- get_g("MR", n)
    mean_C <- get_g("MR_Phage", mean); sd_C <- get_g("MR_Phage", sd); nC <- get_g("MR_Phage", n)
    mean_D <- get_g("MR_AFMT", mean); sd_D <- get_g("MR_AFMT", sd); nD <- get_g("MR_AFMT", n)

    # ANOVA on log2 (if enough data)
    p_overall <- NA_real_
    mse <- NA_real_
    df_res <- NA_real_
    dd <- d %>% filter(!is.na(value_log2))
    if (nrow(dd) >= 8 && length(unique(dd$group)) >= 2) {
      fit <- try(stats::aov(value_log2 ~ group, data = dd), silent = TRUE)
      if (!inherits(fit, "try-error")) {
        s <- summary(fit)[[1]]
        if ("Pr(>F)" %in% colnames(s)) p_overall <- s[["Pr(>F)"]][1]
        # pooled residual variance for LSD
        mse <- s[["Mean Sq"]][nrow(s)]
        df_res <- s[["Df"]][nrow(s)]
      }
    }

    # Pairwise LSD p-values (ANOVA pooled MSE)
    p_A_vs_B <- lsd_p(mean_A, mean_B, nA, nB, mse, df_res)
    p_C_vs_B <- lsd_p(mean_C, mean_B, nC, nB, mse, df_res)
    p_D_vs_B <- lsd_p(mean_D, mean_B, nD, nB, mse, df_res)
    p_C_vs_A <- lsd_p(mean_C, mean_A, nC, nA, mse, df_res)
    p_D_vs_A <- lsd_p(mean_D, mean_A, nD, nA, mse, df_res)
    p_C_vs_D <- lsd_p(mean_C, mean_D, nC, nD, mse, df_res)

    # log2FC
    log2FC_A_vs_B <- mean_A - mean_B
    log2FC_C_vs_B <- mean_C - mean_B
    log2FC_D_vs_B <- mean_D - mean_B
    log2FC_C_vs_A <- mean_C - mean_A
    log2FC_D_vs_A <- mean_D - mean_A
    log2FC_C_vs_D <- mean_C - mean_D

    # geometric means back to original scale
    gm_A <- 2^mean_A
    gm_B <- 2^mean_B
    gm_C <- 2^mean_C
    gm_D <- 2^mean_D

    # rescue relative to GR vs MR baseline, on geometric means
    rescue_phage_pct <- rescue_pct(gm_A, gm_B, gm_C)
    rescue_AFMT_pct  <- rescue_pct(gm_A, gm_B, gm_D)

    tibble::tibble(
      p_overall = p_overall,
      mean_A = mean_A, sd_A = sd_A,
      mean_B = mean_B, sd_B = sd_B,
      mean_C = mean_C, sd_C = sd_C,
      mean_D = mean_D, sd_D = sd_D,

      log2FC_A_vs_B = log2FC_A_vs_B, p_A_vs_B = p_A_vs_B,
      log2FC_C_vs_B = log2FC_C_vs_B, p_C_vs_B = p_C_vs_B,
      log2FC_D_vs_B = log2FC_D_vs_B, p_D_vs_B = p_D_vs_B,
      log2FC_C_vs_A = log2FC_C_vs_A, p_C_vs_A = p_C_vs_A,
      log2FC_D_vs_A = log2FC_D_vs_A, p_D_vs_A = p_D_vs_A,
      log2FC_C_vs_D = log2FC_C_vs_D, p_C_vs_D = p_C_vs_D,

      gm_A = gm_A, gm_B = gm_B, gm_C = gm_C, gm_D = gm_D,
      rescue_phage_pct = rescue_phage_pct,
      rescue_AFMT_pct = rescue_AFMT_pct
    )
  }) %>%
  ungroup()

# Attach meta columns (one row per metabolite)
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
    average_molecular_weight = if ("average Molecular'" %in% names(df_raw)) as_num(`average Molecular'`) else NA_real_
  ) %>%
  distinct()

out <- meta_one %>%
  left_join(stats_tbl, by = "metabolite") %>%
  # BH FDR adjustments
  mutate(
    fdr_overall = p.adjust(p_overall, method = "BH"),
    fdr_A_vs_B = p.adjust(p_A_vs_B, method = "BH"),
    fdr_C_vs_B = p.adjust(p_C_vs_B, method = "BH"),
    fdr_D_vs_B = p.adjust(p_D_vs_B, method = "BH"),
    fdr_C_vs_A = p.adjust(p_C_vs_A, method = "BH"),
    fdr_D_vs_A = p.adjust(p_D_vs_A, method = "BH"),
    fdr_C_vs_D = p.adjust(p_C_vs_D, method = "BH")
  ) %>%
  # Keep a column order similar to platform tables
  select(
    metabolite,
    p_overall,
    mean_A, sd_A, mean_B, sd_B, mean_C, sd_C, mean_D, sd_D,
    log2FC_A_vs_B, p_A_vs_B,
    log2FC_C_vs_B, p_C_vs_B,
    log2FC_D_vs_B, p_D_vs_B,
    log2FC_C_vs_A, p_C_vs_A,
    log2FC_D_vs_A, p_D_vs_A,
    log2FC_C_vs_D, p_C_vs_D,
    fdr_A_vs_B, fdr_C_vs_B, fdr_D_vs_B, fdr_C_vs_A, fdr_D_vs_A, fdr_C_vs_D, fdr_overall,
    CAS, formula, pubchemid, hmdb, kegg, class,
    average_molecular_weight,
    gm_A, gm_B, gm_C, gm_D,
    rescue_phage_pct, rescue_AFMT_pct,
    中文名称, 分类
  )

# ----------------------------
# Write outputs
# ----------------------------
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, OUT_SHEET)
openxlsx::writeData(wb, OUT_SHEET, out)

openxlsx::saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
readr::write_csv(out, sub("\\.xlsx$", ".csv", OUT_XLSX))

message("[DONE] Wrote: ", OUT_XLSX)
message("[DONE] Wrote: ", sub("\\.xlsx$", ".csv", OUT_XLSX))
