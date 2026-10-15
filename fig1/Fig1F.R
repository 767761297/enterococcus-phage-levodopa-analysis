suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(logistf)
})

## =========================================================
## Genus-level unadjusted ROC workflow
## M1: Enterococcus genus abundance
## M2: Enterococcus-host phage burden
## M3: abundance + phage burden
## - No covariate adjustment
## - 16S source: table_filtered_w_tax.tsv.gz
## - Phage definition: all vOTUs with CHERRY host prediction related to Enterococcus
##   (matched by Enterococcus in Host / Host_NCBI_lineage / Host_GTDB_lineage)
## - No CHERRYScore filtering (to stay closest to previous genus-level burden logic)
## =========================================================

## -------------------------
## 0) User settings
## -------------------------
IN_META  <- "metadata.csv"
IN_16S   <- "table_filtered_w_tax.tsv.gz"
IN_VOTU  <- "vOTU_abundance_matrix_TPM.tsv"
IN_CHERY <- "cherry_prediction.tsv"

SEED <- 1
set.seed(SEED)

ROC_COLORS <- c(
  M1_Enterococcus_16S = "#0072B2",
  M2_Enterococcus_phage = "#D55E00",
  M3_Combined = "#009E73"
)

## -------------------------
## 1) Helper functions
## -------------------------
read_metadata_safely <- function(path) {
  tries <- list(
    function() as.data.table(read.csv(path, check.names = FALSE, fileEncoding = "GB18030")),
    function() as.data.table(read.csv(path, check.names = FALSE, fileEncoding = "UTF-8")),
    function() as.data.table(read.csv(path, check.names = FALSE))
  )
  for (fn in tries) {
    out <- tryCatch(fn(), error = function(e) NULL)
    if (!is.null(out) && nrow(out) > 0) return(out)
  }
  stop("无法读取 metadata.csv，请检查文件编码。")
}

roc_points <- function(y, p) {
  y <- as.integer(y); p <- as.numeric(p)
  ok <- is.finite(y) & is.finite(p) & y %in% c(0,1)
  y <- y[ok]; p <- p[ok]
  if (length(unique(y)) < 2) stop("y 只有一个类别，无法计算ROC")

  ord <- order(p, decreasing = TRUE)
  y <- y[ord]; p <- p[ord]

  P <- sum(y == 1); N <- sum(y == 0)
  tp <- cumsum(y == 1)
  fp <- cumsum(y == 0)
  change <- c(TRUE, diff(p) != 0)

  tp_u <- tp[change]
  fp_u <- fp[change]
  thr  <- p[change]

  tpr <- tp_u / P
  fpr <- fp_u / N

  data.table(
    threshold = c(Inf, thr, -Inf),
    FPR = c(0, fpr, 1),
    TPR = c(0, tpr, 1)
  )
}

auc_trapz <- function(roc_dt) {
  d <- roc_dt[order(FPR, TPR)]
  x <- d$FPR; y <- d$TPR
  sum((x[-1] - x[-length(x)]) * (y[-1] + y[-length(y)]) / 2)
}

best_threshold_youden <- function(y, p) {
  r <- roc_points(y, p)
  rr <- r[is.finite(threshold)]
  rr[, youden := TPR - FPR]
  best <- rr[which.max(youden)][1]

  thr <- best$threshold
  pred <- ifelse(p >= thr, 1L, 0L)
  y <- as.integer(y)

  TP <- sum(pred == 1 & y == 1, na.rm = TRUE)
  TN <- sum(pred == 0 & y == 0, na.rm = TRUE)
  FP <- sum(pred == 1 & y == 0, na.rm = TRUE)
  FN <- sum(pred == 0 & y == 1, na.rm = TRUE)

  se  <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  sp  <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  acc <- (TP + TN) / length(y)

  data.table(
    best_threshold = thr,
    sensitivity = se,
    specificity = sp,
    accuracy = acc,
    TP = TP, TN = TN, FP = FP, FN = FN
  )
}

fit_firth <- function(formula, data) {
  ctrl <- logistf.control(maxit = 800, maxstep = 0.25)
  plc  <- logistpl.control(maxit = 800)
  logistf(formula, data = data, control = ctrl, plcontrol = plc)
}

coef_table_from_fit <- function(fit, model_name) {
  b <- coef(fit)
  V <- vcov(fit)
  se <- sqrt(diag(V))

  pv <- suppressWarnings(as.numeric(fit$prob))
  if (length(pv) != length(b) || all(is.na(pv))) {
    z <- as.numeric(b) / as.numeric(se)
    pv <- 2 * pnorm(-abs(z))
  }

  ci <- tryCatch(confint(fit), error = function(e) NULL)
  if (!is.null(ci) && is.matrix(ci) && nrow(ci) == length(b)) {
    lo <- ci[,1]; hi <- ci[,2]
  } else {
    lo <- as.numeric(b) - 1.96 * as.numeric(se)
    hi <- as.numeric(b) + 1.96 * as.numeric(se)
  }

  data.table(
    term = names(b),
    coef = as.numeric(b),
    se = as.numeric(se),
    OR = exp(as.numeric(b)),
    CI_low = exp(as.numeric(lo)),
    CI_high = exp(as.numeric(hi)),
    p = as.numeric(pv),
    model = model_name
  )
}

plot_single_roc <- function(roc_dt, auc_val, model_label, color_hex, out_pdf, out_png) {
  p <- ggplot(roc_dt, aes(FPR, TPR)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey65") +
    geom_step(color = color_hex, linewidth = 1.2) +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
    theme_bw(base_size = 13) +
    labs(
      title = paste0(model_label, " (LOOCV)"),
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    annotate("text", x = 0.98, y = 0.06,
             label = sprintf("AUC = %.4f", auc_val),
             hjust = 1, vjust = 0, size = 4.2)

  ggsave(out_pdf, p, width = 6, height = 5)
  ggsave(out_png, p, width = 6, height = 5, dpi = 300)
}

plot_combined_roc <- function(roc_list, auc_table, out_pdf, out_png) {
  dt <- rbindlist(lapply(names(roc_list), function(nm) {
    x <- copy(roc_list[[nm]])
    x[, model := nm]
    x
  }), fill = TRUE)

  auc_lab <- c(
    sprintf("M1 (AUC=%.4f)", auc_table[model == "M1_Enterococcus_16S", AUC_LOOCV]),
    sprintf("M2 (AUC=%.4f)", auc_table[model == "M2_Enterococcus_phage", AUC_LOOCV]),
    sprintf("M3 (AUC=%.4f)", auc_table[model == "M3_Combined", AUC_LOOCV])
  )
  names(auc_lab) <- c("M1_Enterococcus_16S", "M2_Enterococcus_phage", "M3_Combined")

  p <- ggplot(dt, aes(FPR, TPR, color = model)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey65") +
    geom_step(linewidth = 1.2) +
    scale_color_manual(
      values = ROC_COLORS,
      breaks = c("M1_Enterococcus_16S", "M2_Enterococcus_phage", "M3_Combined"),
      labels = auc_lab
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_bw(base_size = 13) +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      plot.title = element_text(hjust = 0.5)
    ) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
    labs(
      title = "Enterococcus genus ROC (LOOCV)",
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    )

  ggsave(out_pdf, p, width = 8.8, height = 6.2)
  ggsave(out_png, p, width = 8.8, height = 6.2, dpi = 300)
}

## -------------------------
## 2) Output folders
## -------------------------
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR <- paste0("Genus_unadjusted_3models_ROC_", ts)

dir.create(OUT_DIR, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "00_inputs_snapshot"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "01_supporting_tables"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "02_per_sample_table"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "03_model_coefficients"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "04_loocv_predictions"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "05_auc_and_threshold"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "06_roc_points"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "07_plots"), showWarnings = FALSE)

message("输出目录: ", OUT_DIR)

for (ff in c(IN_META, IN_16S, IN_VOTU, IN_CHERY)) {
  if (!file.exists(ff)) stop("缺少输入文件：", ff)
  file.copy(ff, file.path(OUT_DIR, "00_inputs_snapshot", basename(ff)), overwrite = TRUE)
}

## -------------------------
## 3) Read metadata and keep PD_GR / PD_MR
## -------------------------
meta <- read_metadata_safely(IN_META)
need_meta <- c("Sample", "Group", "Levodopa responsiveness (1 MR; 0 GR)")
if (!all(need_meta %in% names(meta))) {
  stop("metadata.csv 缺少关键列：", paste(setdiff(need_meta, names(meta)), collapse = ", "))
}

meta_use <- meta[Group %in% c("PD_GR", "PD_MR"), .(
  Sample = as.character(Sample),
  Group = as.character(Group),
  MR = as.integer(get("Levodopa responsiveness (1 MR; 0 GR)"))
)]
meta_use <- meta_use[MR %in% c(0L, 1L)]
setorder(meta_use, Sample)

if (nrow(meta_use) < 20) stop("PD_GR / PD_MR 有效样本数过少，请检查 metadata.csv")

fwrite(meta_use, file.path(OUT_DIR, "01_supporting_tables", "metadata_PD_GR_MR_used.csv"))

## -------------------------
## 4) Build genus-level Enterococcus 16S counts from table_filtered_w_tax.tsv.gz
## -------------------------
tab16s <- as.data.table(read.delim(gzfile(IN_16S), sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE))
if (!"taxonomy" %in% names(tab16s)) stop("16S 表缺少 taxonomy 列。")

sample_cols_16s <- intersect(meta_use$Sample, names(tab16s))
if (length(sample_cols_16s) == 0) stop("16S 表与 metadata 没有匹配到样本列。")

entero_pat <- "g[_]{1,2}Enterococcus"
asv_keep <- grepl(entero_pat, tab16s$taxonomy, ignore.case = TRUE, perl = TRUE)
entero_asv <- tab16s[asv_keep, c("#ID", sample_cols_16s, "taxonomy"), with = FALSE]

if (nrow(entero_asv) == 0) {
  stop("16S 表中没有匹配到 genus Enterococcus 的 ASV，请检查 taxonomy 格式。")
}

entero16s_counts <- data.table(
  Sample = sample_cols_16s,
  Enterococcus_16S_counts = colSums(as.matrix(entero_asv[, ..sample_cols_16s]), na.rm = TRUE)
)

fwrite(entero_asv, file.path(OUT_DIR, "01_supporting_tables", "Enterococcus_ASVs_from_16S_table_filtered.csv"))
fwrite(entero16s_counts, file.path(OUT_DIR, "01_supporting_tables", "Enterococcus_16S_counts_by_sample.csv"))

## -------------------------
## 5) Build genus-level Enterococcus phage burden from CHERRY + vOTU TPM
## -------------------------
cherry <- fread(IN_CHERY)
votu   <- fread(IN_VOTU)

need_cherry <- c("Accession", "Host", "Host_NCBI_lineage", "Host_GTDB_lineage")
if (!all(need_cherry %in% names(cherry))) {
  stop("cherry_prediction.tsv 缺少关键列：", paste(setdiff(need_cherry, names(cherry)), collapse = ", "))
}
if (!"id" %in% names(votu)) stop("vOTU_abundance_matrix_TPM.tsv 缺少 id 列。")

host_match <- grepl("Enterococcus", cherry$Host, ignore.case = TRUE) |
  grepl("Enterococcus", cherry$Host_NCBI_lineage, ignore.case = TRUE) |
  grepl("Enterococcus", cherry$Host_GTDB_lineage, ignore.case = TRUE)

cherry_entero <- cherry[host_match]
if (nrow(cherry_entero) == 0) stop("CHERRY 结果中没有匹配到 Enterococcus 宿主相关 vOTU。")

entero_ids <- intersect(cherry_entero$Accession, votu$id)
if (length(entero_ids) == 0) stop("CHERRY 的 Enterococcus 相关 Accession 与 vOTU abundance 表没有交集。")

votu_keep <- votu[id %in% entero_ids]
sample_cols_votu <- intersect(meta_use$Sample, names(votu_keep))
if (length(sample_cols_votu) == 0) stop("vOTU abundance 表与 metadata 没有匹配到样本列。")

entero_phage_tpm <- data.table(
  Sample = sample_cols_votu,
  Enterococcus_phage_burden_TPM = colSums(as.matrix(votu_keep[, ..sample_cols_votu]), na.rm = TRUE)
)

fwrite(cherry_entero, file.path(OUT_DIR, "01_supporting_tables", "Enterococcus_host_related_vOTUs_from_CHERRY.csv"))
fwrite(votu_keep, file.path(OUT_DIR, "01_supporting_tables", "Enterococcus_host_related_vOTU_TPM_matrix.csv"))
fwrite(entero_phage_tpm, file.path(OUT_DIR, "01_supporting_tables", "Enterococcus_phage_burden_TPM_by_sample.csv"))

## -------------------------
## 6) Merge per-sample table
## -------------------------
ps <- merge(meta_use, entero16s_counts, by = "Sample", all.x = TRUE)
ps <- merge(ps, entero_phage_tpm, by = "Sample", all.x = TRUE)

ps[is.na(Enterococcus_16S_counts), Enterococcus_16S_counts := 0]
ps[is.na(Enterococcus_phage_burden_TPM), Enterococcus_phage_burden_TPM := 0]
ps[, log_Enterococcus_16S := log1p(Enterococcus_16S_counts)]
ps[, log_Enterococcus_phage := log1p(Enterococcus_phage_burden_TPM)]

setcolorder(ps, c(
  "Sample", "Group", "MR",
  "Enterococcus_16S_counts", "log_Enterococcus_16S",
  "Enterococcus_phage_burden_TPM", "log_Enterococcus_phage"
))

fwrite(ps, file.path(OUT_DIR, "02_per_sample_table", "per_sample_Genus_unadjusted_16S_plus_phage_table.csv"))

## -------------------------
## 7) Full-sample fit (for coefficient tables)
## -------------------------
db <- copy(ps)
db[, z_abund := if (sd(log_Enterococcus_16S) == 0) 0 else as.numeric(scale(log_Enterococcus_16S))]
db[, z_phage := if (sd(log_Enterococcus_phage) == 0) 0 else as.numeric(scale(log_Enterococcus_phage))]

models <- list(
  M1_Enterococcus_16S   = MR ~ z_abund,
  M2_Enterococcus_phage = MR ~ z_phage,
  M3_Combined           = MR ~ z_abund + z_phage
)

coef_all <- list()
for (nm in names(models)) {
  fit <- fit_firth(models[[nm]], db)

  sink(file.path(OUT_DIR, "03_model_coefficients", paste0(nm, "_firth_summary.txt")))
  print(summary(fit))
  sink()

  ct <- coef_table_from_fit(fit, nm)
  fwrite(ct, file.path(OUT_DIR, "03_model_coefficients", paste0(nm, "_coefficients.csv")))
  coef_all[[nm]] <- ct
}
fwrite(rbindlist(coef_all, fill = TRUE),
       file.path(OUT_DIR, "03_model_coefficients", "Genus_unadjusted_Firth_coefficients_all_models.csv"))

## -------------------------
## 8) LOOCV prediction (train-fold standardization only)
## -------------------------
make_loocv_pred <- function(model_key, data) {
  n <- nrow(data)
  pred <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    tr <- copy(data[-i])
    te <- copy(data[i])

    m_ab <- mean(tr$log_Enterococcus_16S); s_ab <- sd(tr$log_Enterococcus_16S)
    if (!is.finite(s_ab) || s_ab == 0) s_ab <- 1
    m_ph <- mean(tr$log_Enterococcus_phage); s_ph <- sd(tr$log_Enterococcus_phage)
    if (!is.finite(s_ph) || s_ph == 0) s_ph <- 1

    tr[, z_abund := (log_Enterococcus_16S - m_ab) / s_ab]
    te[, z_abund := (log_Enterococcus_16S - m_ab) / s_ab]
    tr[, z_phage := (log_Enterococcus_phage - m_ph) / s_ph]
    te[, z_phage := (log_Enterococcus_phage - m_ph) / s_ph]

    fml <- switch(
      model_key,
      "M1_Enterococcus_16S"   = MR ~ z_abund,
      "M2_Enterococcus_phage" = MR ~ z_phage,
      "M3_Combined"           = MR ~ z_abund + z_phage
    )

    fit <- tryCatch(fit_firth(fml, tr), error = function(e) NULL)
    if (is.null(fit)) next

    pred[i] <- tryCatch(as.numeric(predict(fit, newdata = te, type = "response")),
                        error = function(e) NA_real_)
  }
  pred
}

pred_tab <- ps[, .(
  Sample, Group, MR,
  Enterococcus_16S_counts, log_Enterococcus_16S,
  Enterococcus_phage_burden_TPM, log_Enterococcus_phage
)]

for (nm in names(models)) {
  message("LOOCV: ", nm)
  pred_tab[[paste0("pred_", nm)]] <- make_loocv_pred(nm, pred_tab)
}

fwrite(pred_tab, file.path(OUT_DIR, "04_loocv_predictions", "Genus_unadjusted_LOOCV_per_sample_predictions.csv"))

## -------------------------
## 9) ROC / AUC / thresholds
## -------------------------
pretty_map <- c(
  M1_Enterococcus_16S = "M1: Enterococcus genus abundance",
  M2_Enterococcus_phage = "M2: Enterococcus genus phage burden",
  M3_Combined = "M3: abundance + phage burden"
)

roc_list <- list()
auc_rows <- list()
threshold_rows <- list()

for (nm in names(models)) {
  pc <- paste0("pred_", nm)
  rp <- roc_points(pred_tab$MR, pred_tab[[pc]])
  auc_val <- auc_trapz(rp)
  bt <- best_threshold_youden(pred_tab$MR, pred_tab[[pc]])

  roc_list[[nm]] <- rp
  auc_rows[[nm]] <- data.table(
    model = nm,
    pretty_label = pretty_map[[nm]],
    AUC_LOOCV = auc_val
  )
  threshold_rows[[nm]] <- cbind(
    data.table(model = nm, pretty_label = pretty_map[[nm]], AUC_LOOCV = auc_val),
    bt
  )

  rp_out <- copy(rp)
  rp_out[, model := nm]
  rp_out[, pretty_label := pretty_map[[nm]]]
  fwrite(rp_out, file.path(OUT_DIR, "06_roc_points", paste0("roc_points_", nm, ".csv")))
}

auc_tab <- rbindlist(auc_rows)
metric_tab <- rbindlist(threshold_rows)

# choose best model: highest AUC, then highest accuracy, then highest sensitivity+specificity
metric_tab[, tie_key := sensitivity + specificity]
setorder(metric_tab, -AUC_LOOCV, -accuracy, -tie_key)
metric_tab[, rank_by_diagnostic_performance := seq_len(.N)]
metric_tab[, best_model_flag := ifelse(rank_by_diagnostic_performance == 1L, "YES", "NO")]

best_model <- metric_tab[1, model]
best_pretty <- metric_tab[1, pretty_label]

fwrite(metric_tab, file.path(OUT_DIR, "05_auc_and_threshold", "Genus_unadjusted_LOOCV_AUC_threshold_summary.csv"))
fwrite(metric_tab[best_model_flag == "YES"],
       file.path(OUT_DIR, "05_auc_and_threshold", "Best_model_selected_for_followup.csv"))

## -------------------------
## 10) ROC plots
## -------------------------
for (nm in names(models)) {
  plot_single_roc(
    roc_dt = roc_list[[nm]],
    auc_val = auc_tab[model == nm, AUC_LOOCV][1],
    model_label = pretty_map[[nm]],
    color_hex = ROC_COLORS[[nm]],
    out_pdf = file.path(OUT_DIR, "07_plots", paste0("ROC_", nm, ".pdf")),
    out_png = file.path(OUT_DIR, "07_plots", paste0("ROC_", nm, ".png"))
  )
}

# combined plot order = M1, M2, M3
auc_plot <- auc_tab[match(names(models), model)]
plot_combined_roc(
  roc_list = roc_list,
  auc_table = auc_plot,
  out_pdf = file.path(OUT_DIR, "07_plots", "ROC_M1_M2_M3_combined.pdf"),
  out_png = file.path(OUT_DIR, "07_plots", "ROC_M1_M2_M3_combined.png")
)

## -------------------------
## 11) concise readme-style summary
## -------------------------
summary_txt <- c(
  paste0("Output directory: ", OUT_DIR),
  "",
  "Core analysis design:",
  "- No covariate adjustment",
  "- M1 = Enterococcus genus abundance only",
  "- M2 = Enterococcus genus phage burden only",
  "- M3 = abundance + phage burden",
  "- LOOCV-based predicted probabilities used for ROC evaluation",
  "",
  "16S definition:",
  "- table_filtered_w_tax.tsv.gz",
  "- all ASVs with taxonomy containing g_Enterococcus aggregated to sample-level counts",
  "",
  "Phage definition:",
  "- all CHERRY-predicted Enterococcus-related host vOTUs (matched by Enterococcus in Host / lineage fields)",
  "- sample-level burden = sum of TPM across matched vOTUs",
  "- no CHERRYScore filter applied",
  "",
  paste0("Best model selected for follow-up: ", best_model, "  [", best_pretty, "]"),
  paste0("Best model AUC_LOOCV: ", sprintf("%.4f", metric_tab[1, AUC_LOOCV])),
  paste0("Best threshold: ", sprintf("%.6f", metric_tab[1, best_threshold])),
  paste0("Sensitivity: ", sprintf("%.4f", metric_tab[1, sensitivity])),
  paste0("Specificity: ", sprintf("%.4f", metric_tab[1, specificity])),
  paste0("Accuracy: ", sprintf("%.4f", metric_tab[1, accuracy]))
)
writeLines(summary_txt, con = file.path(OUT_DIR, "README_summary.txt"))

message("\n✅ 完成：unadjusted 三模型 ROC 脚本已运行到输出阶段。")
message("优先查看：")
message("  1) 02_per_sample_table/per_sample_Genus_unadjusted_16S_plus_phage_table.csv")
message("  2) 05_auc_and_threshold/Genus_unadjusted_LOOCV_AUC_threshold_summary.csv")
message("  3) 05_auc_and_threshold/Best_model_selected_for_followup.csv")
message("  4) 07_plots/ROC_M1_M2_M3_combined.png")
message("  5) 07_plots/ROC_M1_Enterococcus_16S.png")
message("  6) 07_plots/ROC_M2_Enterococcus_phage.png")
message("  7) 07_plots/ROC_M3_Combined.png")