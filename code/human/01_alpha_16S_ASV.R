## =========================
## 01_alpha_16S_ASV.R
## =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(readr)
  library(janitor)
  library(vegan)
  library(ggplot2)
})

## ---------- helpers ----------
detect_file <- function(prefer, patterns) {
  if (!is.null(prefer) && file.exists(prefer)) return(prefer)
  files <- list.files(".", full.names = FALSE)
  for (pat in patterns) {
    hit <- files[str_detect(files, regex(pat, ignore_case = TRUE))]
    if (length(hit) > 0) return(hit[1])
  }
  stop("找不到输入文件。请检查工作目录文件名是否匹配：", paste(patterns, collapse = " | "))
}

read_csv_smart <- function(path) {
  tryCatch(
    readr::read_csv(path, show_col_types = FALSE, locale = locale(encoding = "GB18030")),
    error = function(e) readr::read_csv(path, show_col_types = FALSE, locale = locale(encoding = "UTF-8"))
  )
}

normalize_id <- function(x) {
  x |> as.character() |> str_trim() |> str_replace_all("\\s+", "")
}

find_col <- function(df, patterns) {
  nms <- names(df)
  for (p in patterns) {
    hit <- nms[str_detect(nms, regex(p, ignore_case = TRUE))]
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

write_tsv0 <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(x, path)
}

## ---------- output dirs ----------
out_root <- "results_16S_ASV_alpha"
dir.create(out_root, showWarnings = FALSE)
dir.create(file.path(out_root, "tables"), showWarnings = FALSE)
dir.create(file.path(out_root, "figures"), showWarnings = FALSE)
dir.create(file.path(out_root, "logs"), showWarnings = FALSE)
dir.create(file.path(out_root, "00_manifest"), showWarnings = FALSE)

## ---------- locate inputs (no paths) ----------
meta_file <- detect_file(
  prefer   = "metadata.csv",
  patterns = c("^metadata\\.csv$")
)

asv_file <- detect_file(
  prefer   = "table_f_even_depth_w_tax.xls",
  patterns = c("table_f_even_depth_w_tax\\.(xls|tsv|txt)$")
)

## ---------- load metadata ----------
meta_raw <- read_csv_smart(meta_file) |> janitor::clean_names()

sample_col <- find_col(meta_raw, c("^sample$"))
group_col  <- find_col(meta_raw, c("^group$"))
age_col    <- find_col(meta_raw, c("^age$"))
bmi_col    <- find_col(meta_raw, c("^bmi"))
dur_col    <- find_col(meta_raw, c("disease_duration"))
con_col    <- find_col(meta_raw, c("constipation"))
ldopa_col  <- find_col(meta_raw, c("levodopa_daily_dose"))

stopifnot(!is.na(sample_col), !is.na(group_col), !is.na(age_col), !is.na(bmi_col),
          !is.na(dur_col), !is.na(con_col))

meta <- meta_raw |>
  mutate(sample_id = normalize_id(.data[[sample_col]]),
         group = .data[[group_col]]) |>
  filter(!is.na(sample_id) & sample_id != "") |>                # 去掉末尾空白行
  filter(group %in% c("PD_GR","PD_MR")) |>                       # 仅MR/GR
  transmute(
    sample_id,
    group = factor(group, levels = c("PD_GR","PD_MR")),
    age = as.numeric(.data[[age_col]]),
    bmi = as.numeric(.data[[bmi_col]]),
    disease_duration = as.numeric(.data[[dur_col]]),
    constipation = factor(as.integer(.data[[con_col]]), levels = c(0,1), labels = c("No","Yes")),
    levodopa_dose = if (!is.na(ldopa_col)) as.numeric(.data[[ldopa_col]]) else NA_real_
  )

stopifnot(nrow(meta) == 40)

write_tsv0(meta, file.path(out_root, "tables", "metadata_model_40_PD_MR_GR.tsv"))

## ---------- load 16S ASV table ----------
## 注意：该“xls”常为tab文本；优先用fread读取
asv_dt <- data.table::fread(asv_file, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)

## manifest
manifest <- tibble(
  file = c(meta_file, asv_file),
  rows = c(nrow(meta_raw), nrow(asv_dt)),
  cols = c(ncol(meta_raw), ncol(asv_dt))
)
write_tsv0(manifest, file.path(out_root, "00_manifest", "inputs_manifest.tsv"))

## identify columns
id_col  <- if ("#ID" %in% names(asv_dt)) "#ID" else names(asv_dt)[1]
tax_col <- if ("taxonomy" %in% names(asv_dt)) "taxonomy" else NA_character_

sample_cols <- names(asv_dt)[str_detect(names(asv_dt), "^PD\\d+$")]
stopifnot(length(sample_cols) > 0)

## keep only the 40 samples (intersection; order = metadata)
keep_samples <- intersect(meta$sample_id, sample_cols)
if (length(keep_samples) != 40) {
  stop("16S表与metadata样本交集不是40。请检查样本列名是否为PDxx并与metadata$Sample一致。")
}
keep_samples <- meta$sample_id  # enforce order

## build matrix: features x samples
mat <- asv_dt[, keep_samples, drop = FALSE]
mat[] <- lapply(mat, as.numeric)
mat <- as.matrix(mat)

## basic QC tables (per-sample)
lib_size <- colSums(mat, na.rm = TRUE)
observed <- colSums(mat > 0, na.rm = TRUE)

qc_per_sample <- meta |>
  mutate(
    library_size = lib_size[sample_id],
    observed_asv = observed[sample_id]
  )
write_tsv0(qc_per_sample, file.path(out_root, "tables", "QC_per_sample_librarysize_observed.tsv"))

## alpha diversity metrics
comm <- t(mat)  # samples x features (40 x ~200k)
shannon <- vegan::diversity(comm, index = "shannon")
simpson <- vegan::diversity(comm, index = "simpson")  # 1-lambda

alpha_wide <- meta |>
  mutate(
    observed = observed[sample_id],
    shannon  = as.numeric(shannon[sample_id]),
    simpson  = as.numeric(simpson[sample_id]),
    inv_simpson = 1 / (1 - simpson),
    exp_shannon = exp(shannon),
    pielou_evenness = ifelse(observed > 1, shannon / log(observed), NA_real_)
# optional; interpret carefully
  )

write_tsv0(alpha_wide, file.path(out_root, "tables", "alpha_16S_per_sample_wide.tsv"))

alpha_long <- alpha_wide |>
  select(sample_id, group, age, bmi, disease_duration, constipation, levodopa_dose,
         observed, shannon, simpson, inv_simpson, exp_shannon, pielou_evenness) |>
  pivot_longer(cols = c(observed, shannon, simpson, inv_simpson, exp_shannon, pielou_evenness),
               names_to = "metric", values_to = "value")
write_tsv0(alpha_long, file.path(out_root, "tables", "alpha_16S_per_sample_long.tsv"))

## ---------- stats: unadjusted + primary adjusted + sensitivity(add dose) ----------
wilcox_one <- function(df, metric) {
  x <- df |> filter(metric == !!metric)
  w <- wilcox.test(value ~ group, data = x, exact = FALSE)
  tibble(metric = metric, test = "wilcoxon", p_value = w$p.value)
}

lm_one <- function(df, metric, add_dose = FALSE) {
  x <- df |> filter(metric == !!metric)
  if (add_dose && all(is.na(x$levodopa_dose))) return(NULL)

  fml <- if (!add_dose) {
    value ~ group + age + bmi + disease_duration + constipation
  } else {
    value ~ group + age + bmi + disease_duration + constipation + levodopa_dose
  }
  fit <- lm(fml, data = x)
  s <- summary(fit)$coefficients
  ## groupPD_MR effect (PD_MR vs PD_GR)
  rn <- rownames(s)
  grp_row <- rn[str_detect(rn, "^groupPD_MR$")]
  if (length(grp_row) == 0) grp_row <- rn[str_detect(rn, "^group")]
  eff <- s[grp_row[1], , drop = FALSE]

  tibble(
    metric = metric,
    model = ifelse(add_dose, "adjusted_sensitivity_add_dose", "adjusted_primary"),
    beta_group_PD_MR = eff[1, "Estimate"],
    se = eff[1, "Std. Error"],
    t = eff[1, "t value"],
    p_value = eff[1, "Pr(>|t|)"]
  )
}

metrics <- unique(alpha_long$metric)

stats_unadj <- bind_rows(lapply(metrics, \(m) wilcox_one(alpha_long, m)))
stats_primary <- bind_rows(lapply(metrics, \(m) lm_one(alpha_long, m, add_dose = FALSE)))
stats_sens <- bind_rows(lapply(metrics, \(m) lm_one(alpha_long, m, add_dose = TRUE)))

write_tsv0(stats_unadj,  file.path(out_root, "tables", "alpha_stats_unadjusted_wilcoxon.tsv"))
write_tsv0(stats_primary, file.path(out_root, "tables", "alpha_stats_adjusted_primary.tsv"))
if (!is.null(stats_sens) && nrow(stats_sens) > 0) {
  write_tsv0(stats_sens, file.path(out_root, "tables", "alpha_stats_adjusted_sensitivity_add_dose.tsv"))
}

## ---------- minimal figures (per-metric boxplot + points) ----------
pdf(file.path(out_root, "figures", "alpha_16S_boxplots.pdf"), width = 7, height = 5)
for (m in metrics) {
  p <- alpha_long |>
    filter(metric == m) |>
    ggplot(aes(x = group, y = value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, height = 0, size = 2) +
    theme_bw() +
    labs(title = paste0("16S alpha: ", m), x = NULL, y = m)
  print(p)
}
dev.off()

## ---------- logs ----------
capture.output(sessionInfo(), file = file.path(out_root, "logs", "sessionInfo.txt"))
message("16S ASV alpha done -> ", out_root)