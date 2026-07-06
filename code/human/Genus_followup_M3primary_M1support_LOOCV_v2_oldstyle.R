suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(logistf)
})

## =========================================================
## Genus-level follow-up model evaluation (LOOCV-based)
## Primary model:   M3 = Enterococcus genus abundance + phage burden
## Support model:   M1 = Enterococcus genus abundance only
## - No covariate adjustment
## - ROC / Calibration / DCA / Nomogram
## - Single-model figures for M3 and M1
## - Combined comparison figures for M3 vs M1
## - 16S source: table_filtered_w_tax.xls
## - Phage definition: all CHERRY-predicted Enterococcus-related host vOTUs
## =========================================================

## -------------------------
## 0) User settings
## -------------------------
IN_META  <- "metadata.csv"
IN_16S   <- "table_filtered_w_tax.xls"
IN_VOTU  <- "vOTU_abundance_matrix_TPM.tsv"
IN_CHERY <- "cherry_prediction.tsv"

SEED <- 1
set.seed(SEED)

NBINS_CAL <- 8
EPS <- 1e-6
DCA_THRESHOLDS <- seq(0.01, 0.99, by = 0.01)

MODEL_COLORS <- c(
  M1_Enterococcus_16S = "#0072B2",
  M3_Combined = "#009E73",
  Treat_all = "#D55E00",
  Treat_none = "#7F7F7F"
)

PRETTY_MAP <- c(
  M1_Enterococcus_16S = "M1: Enterococcus genus abundance",
  M3_Combined = "M3: abundance + phage burden"
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

clamp_prob <- function(p, eps = EPS) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}

roc_points <- function(y, p) {
  y <- as.integer(y); p <- as.numeric(p)
  ok <- is.finite(y) & is.finite(p) & y %in% c(0, 1)
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
  pred <- ifelse(as.numeric(p) >= thr, 1L, 0L)
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

fit_firth_stable <- function(formula, data) {
  ctrl <- logistf.control(maxit = 5000, maxstep = 0.5)
  logistf(formula, data = data, firth = TRUE, pl = FALSE, control = ctrl)
}

coef_table_from_fit <- function(fit, model_name) {
  b <- coef(fit)
  V <- vcov(fit)
  se <- sqrt(diag(V))
  z  <- as.numeric(b) / as.numeric(se)
  pv <- 2 * pnorm(-abs(z))

  out <- data.table(
    term = names(b),
    coef = as.numeric(b),
    se = as.numeric(se),
    OR = exp(as.numeric(b)),
    CI_low = exp(as.numeric(b) - 1.96 * as.numeric(se)),
    CI_high = exp(as.numeric(b) + 1.96 * as.numeric(se)),
    p_value_wald = as.numeric(pv),
    model = model_name
  )
  out[]
}

safe_train_scale <- function(tr, te, col_raw, col_z) {
  m <- mean(tr[[col_raw]])
  s <- sd(tr[[col_raw]])
  if (!is.finite(s) || s == 0) s <- 1
  tr[[col_z]] <- (tr[[col_raw]] - m) / s
  te[[col_z]] <- (te[[col_raw]] - m) / s
  list(tr = tr, te = te, mean = m, sd = s)
}

make_loocv_pred <- function(model_key, data) {
  n <- nrow(data)
  pred <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    tr <- copy(data[-i])
    te <- copy(data[i])

    tmp <- safe_train_scale(tr, te, "log_Enterococcus_16S", "z_abund")
    tr <- tmp$tr; te <- tmp$te

    if (model_key == "M3_Combined") {
      tmp <- safe_train_scale(tr, te, "log_Enterococcus_phage", "z_phage")
      tr <- tmp$tr; te <- tmp$te
    }

    fml <- switch(
      model_key,
      "M1_Enterococcus_16S" = MR ~ z_abund,
      "M3_Combined" = MR ~ z_abund + z_phage
    )

    fit <- tryCatch(fit_firth_stable(fml, tr), error = function(e) NULL)
    if (is.null(fit)) next

    pred[i] <- tryCatch(as.numeric(predict(fit, newdata = te, type = "response")),
                        error = function(e) NA_real_)
  }

  clamp_prob(pred)
}

make_calibration_table <- function(y, p, model_name, nbins = NBINS_CAL) {
  dt <- data.table(MR = as.integer(y), p = clamp_prob(p))
  dt <- dt[!is.na(MR) & !is.na(p)]
  dt[, r := frank(p, ties.method = "average")]
  dt[, bin := pmax(1L, pmin(nbins, ceiling(r / (.N / nbins))))]

  cal_tab <- dt[, .(
    n = .N,
    mean_pred = mean(p),
    obs_rate = mean(MR)
  ), by = bin][order(mean_pred)]

  cal_tab[, se_obs := sqrt(pmax(obs_rate * (1 - obs_rate) / n, 0))]
  cal_tab[, ci_low := pmax(0, obs_rate - 1.96 * se_obs)]
  cal_tab[, ci_high := pmin(1, obs_rate + 1.96 * se_obs)]
  cal_tab[, model := model_name]
  cal_tab[]
}

calibration_stats <- function(y, p, model_name) {
  dt <- data.table(MR = as.integer(y), p = clamp_prob(p))
  dt[, lp := log(p / (1 - p))]
  fit <- glm(MR ~ lp, data = dt, family = binomial)
  data.table(
    model = model_name,
    intercept = unname(coef(fit)[1]),
    slope = unname(coef(fit)[2]),
    brier = mean((dt$p - dt$MR)^2),
    n = nrow(dt),
    nbins = NBINS_CAL
  )
}

make_dca_table <- function(y, p, model_name, thresholds = DCA_THRESHOLDS) {
  dt <- data.table(MR = as.integer(y), p = clamp_prob(p))
  N <- nrow(dt)
  prev <- mean(dt$MR)

  out <- rbindlist(lapply(thresholds, function(t) {
    pred_pos <- dt$p >= t
    TP <- sum(pred_pos & dt$MR == 1)
    FP <- sum(pred_pos & dt$MR == 0)

    NB_model <- (TP / N) - (FP / N) * (t / (1 - t))
    NB_all   <- prev - (1 - prev) * (t / (1 - t))
    NB_none  <- 0

    data.table(
      threshold = t,
      NB_model = NB_model,
      NB_treat_all = NB_all,
      NB_treat_none = NB_none,
      model = model_name
    )
  }))
  out[]
}

plot_roc_single <- function(roc_dt, auc_val, title_txt, color_hex, out_pdf, out_png) {
  p <- ggplot(roc_dt, aes(FPR, TPR)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey65") +
    geom_step(color = color_hex, linewidth = 1.2) +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
    theme_bw(base_size = 13) +
    labs(
      title = title_txt,
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    annotate("text", x = 0.98, y = 0.06,
             label = sprintf("AUC = %.4f", auc_val),
             hjust = 1, vjust = 0, size = 4.2)

  ggsave(out_pdf, p, width = 6.0, height = 5.0)
  ggsave(out_png, p, width = 6.0, height = 5.0, dpi = 300)
}

plot_roc_combined <- function(roc_m1, roc_m3, auc_m1, auc_m3, out_pdf, out_png) {
  dt <- rbindlist(list(
    copy(roc_m1)[, model := "M1_Enterococcus_16S"],
    copy(roc_m3)[, model := "M3_Combined"]
  ))

  lab <- c(
    sprintf("M1 (AUC=%.4f)", auc_m1),
    sprintf("M3 (AUC=%.4f)", auc_m3)
  )
  names(lab) <- c("M1_Enterococcus_16S", "M3_Combined")

  p <- ggplot(dt, aes(FPR, TPR, color = model)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey65") +
    geom_step(linewidth = 1.2) +
    scale_color_manual(
      values = MODEL_COLORS[c("M1_Enterococcus_16S", "M3_Combined")],
      breaks = c("M1_Enterococcus_16S", "M3_Combined"),
      labels = lab
    ) +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
    theme_bw(base_size = 13) +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      plot.title = element_text(hjust = 0.5)
    ) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
    labs(
      title = "Enterococcus genus ROC (LOOCV): M3 vs M1",
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    )

  ggsave(out_pdf, p, width = 8.2, height = 6.1)
  ggsave(out_png, p, width = 8.2, height = 6.1, dpi = 300)
}

plot_calibration_single <- function(cal_tab, stats_dt, title_txt, color_hex, out_pdf, out_png) {
  p <- ggplot(cal_tab, aes(mean_pred, obs_rate)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey65") +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0, color = "black", linewidth = 0.6) +
    geom_line(color = "black", linewidth = 0.8) +
    geom_point(aes(size = n), color = "black") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    ) +
    labs(
      title = title_txt,
      x = "Predicted probability (LOOCV)",
      y = "Observed MR rate",
      size = "Bin n"
    )

  ggsave(out_pdf, p, width = 6.2, height = 5.2)
  ggsave(out_png, p, width = 6.2, height = 5.2, dpi = 300)
}

plot_calibration_combined <- function(cal_m1, cal_m3, out_pdf, out_png) {
  dt <- rbindlist(list(cal_m3, cal_m1), fill = TRUE)
  dt[, model := factor(model, levels = c("M3_Combined", "M1_Enterococcus_16S"))]

  p <- ggplot(dt, aes(mean_pred, obs_rate, color = model, group = model)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey65") +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0, linewidth = 0.6, alpha = 0.75) +
    geom_line(linewidth = 0.9) +
    geom_point(aes(size = n)) +
    scale_color_manual(
      values = c(
        M3_Combined = "black",
        M1_Enterococcus_16S = "#0072B2"
      ),
      breaks = c("M3_Combined", "M1_Enterococcus_16S"),
      labels = c("M3 model", "M1 model")
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_bw(base_size = 14) +
    theme(
      legend.title = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 14)
    ) +
    labs(
      title = "Calibration curve (LOOCV-based): M3 vs M1",
      x = "Predicted probability (LOOCV)",
      y = "Observed MR rate",
      size = "Bin n"
    )

  ggsave(out_pdf, p, width = 7.2, height = 5.6)
  ggsave(out_png, p, width = 7.2, height = 5.6, dpi = 300)
}

plot_dca_single <- function(dca_tab, title_txt, model_color, out_pdf, out_png) {
  dt <- melt(
    dca_tab[, .(threshold, NB_model, NB_treat_all, NB_treat_none)],
    id.vars = "threshold",
    variable.name = "strategy",
    value.name = "net_benefit"
  )

  p <- ggplot(dt, aes(threshold, net_benefit, color = strategy)) +
    geom_line(linewidth = 1.0) +
    scale_color_manual(
      values = c(
        NB_model = "black",
        NB_treat_all = "#D55E00",
        NB_treat_none = "#0072B2"
      ),
      breaks = c("NB_model", "NB_treat_all", "NB_treat_none"),
      labels = c(
        NB_model = "Model",
        NB_treat_all = "Treat all",
        NB_treat_none = "Treat none"
      )
    ) +
    coord_cartesian(ylim = c(-0.05, 0.5)) +
    theme_bw(base_size = 14) +
    theme(
      legend.title = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 14)
    ) +
    labs(
      title = title_txt,
      x = "Threshold probability",
      y = "Net benefit"
    )

  ggsave(out_pdf, p, width = 7.0, height = 5.8)
  ggsave(out_png, p, width = 7.0, height = 5.8, dpi = 300)
}

plot_dca_combined <- function(dca_m1, dca_m3, out_pdf, out_png) {
  dt_m1 <- copy(dca_m1)[, .(threshold, net_benefit = NB_model, strategy = "M1_Enterococcus_16S")]
  dt_m3 <- copy(dca_m3)[, .(threshold, net_benefit = NB_model, strategy = "M3_Combined")]
  dt_ref <- copy(dca_m3)[, .(
    threshold,
    NB_treat_all,
    NB_treat_none
  )]
  dt_ref <- melt(dt_ref, id.vars = "threshold", variable.name = "strategy", value.name = "net_benefit")

  dt <- rbindlist(list(dt_m3, dt_m1, dt_ref), fill = TRUE)

  p <- ggplot(dt, aes(threshold, net_benefit, color = strategy)) +
    geom_line(linewidth = 1.0) +
    scale_color_manual(
      values = c(
        M3_Combined = "black",
        M1_Enterococcus_16S = "#009E73",
        NB_treat_all = "#D55E00",
        NB_treat_none = "#0072B2"
      ),
      breaks = c("M3_Combined", "M1_Enterococcus_16S", "NB_treat_all", "NB_treat_none"),
      labels = c(
        M3_Combined = "M3 model",
        M1_Enterococcus_16S = "M1 model",
        NB_treat_all = "Treat all",
        NB_treat_none = "Treat none"
      )
    ) +
    coord_cartesian(ylim = c(-0.05, 0.5)) +
    theme_bw(base_size = 14) +
    theme(
      legend.title = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 14)
    ) +
    labs(
      title = "Decision Curve Analysis (LOOCV-based): M3 vs M1",
      x = "Threshold probability",
      y = "Net benefit"
    )

  ggsave(out_pdf, p, width = 7.4, height = 5.8)
  ggsave(out_png, p, width = 7.4, height = 5.8, dpi = 300)
}

make_nomogram_axes <- function(fit, dat, model_name) {
  beta <- coef(fit)
  b0 <- unname(beta["(Intercept)"])
  b1 <- unname(beta["z_abund"])
  if (!is.finite(b1)) stop("Nomogram: z_abund 系数不存在。")

  z_ab_range <- range(dat$z_abund, na.rm = TRUE)

  if (model_name == "M3_Combined") {
    b2 <- unname(beta["z_phage"])
    if (!is.finite(b2)) stop("Nomogram: z_phage 系数不存在。")
    z_ph_range <- range(dat$z_phage, na.rm = TRUE)
    contrib <- c(b1 * z_ab_range, b2 * z_ph_range)
  } else {
    b2 <- NA_real_
    z_ph_range <- c(NA_real_, NA_real_)
    contrib <- c(b1 * z_ab_range)
  }

  max_abs_contrib <- max(abs(contrib), na.rm = TRUE)
  if (!is.finite(max_abs_contrib) || max_abs_contrib == 0) max_abs_contrib <- 1
  scale_points <- 100 / max_abs_contrib

  ax_ab <- data.table(
    var = "z_abund",
    z = seq(z_ab_range[1], z_ab_range[2], length.out = 200),
    points = b1 * seq(z_ab_range[1], z_ab_range[2], length.out = 200) * scale_points
  )

  if (model_name == "M3_Combined") {
    ax_ph <- data.table(
      var = "z_phage",
      z = seq(z_ph_range[1], z_ph_range[2], length.out = 200),
      points = b2 * seq(z_ph_range[1], z_ph_range[2], length.out = 200) * scale_points
    )
    tp <- seq(min(ax_ab$points) + min(ax_ph$points),
              max(ax_ab$points) + max(ax_ph$points), length.out = 300)
  } else {
    ax_ph <- NULL
    tp <- seq(min(ax_ab$points), max(ax_ab$points), length.out = 300)
  }

  ax_total <- data.table(
    total_points = tp,
    p = plogis(b0 + tp / scale_points)
  )

  list(
    ax_ab = ax_ab,
    ax_ph = ax_ph,
    ax_total = ax_total,
    scale_points = scale_points,
    intercept = b0,
    beta_ab = b1,
    beta_ph = b2
  )
}

plot_nomogram_axis <- function(ax_dt, title_txt, xlab_txt, out_pdf, out_png) {
  p <- ggplot(ax_dt, aes(z, points)) +
    geom_line(linewidth = 1.1) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(title = title_txt, x = xlab_txt, y = "Points")

  ggsave(out_pdf, p, width = 6.5, height = 4.6)
  ggsave(out_png, p, width = 6.5, height = 4.6, dpi = 300)
}

plot_nomogram_total <- function(ax_dt, title_txt, out_pdf, out_png) {
  p <- ggplot(ax_dt, aes(total_points, p)) +
    geom_line(linewidth = 1.1) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(title = title_txt, x = "Total points", y = "Predicted probability")

  ggsave(out_pdf, p, width = 6.5, height = 4.6)
  ggsave(out_png, p, width = 6.5, height = 4.6, dpi = 300)
}

## -------------------------
## 2) Output folders
## -------------------------
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR <- paste0("Genus_followup_M3primary_M1support_", ts)

dir.create(OUT_DIR, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "00_inputs_snapshot"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "01_supporting_tables"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "02_per_sample_table"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "03_model_coefficients"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "04_loocv_predictions"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "05_metrics"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "06_roc_points"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "07_calibration_tables"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "08_dca_tables"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "09_nomogram_points"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "10_plots_M3_primary"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "11_plots_M1_support"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "12_plots_M3_vs_M1"), showWarnings = FALSE)

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
## 4) Build genus-level Enterococcus 16S counts from table_filtered_w_tax.xls
## -------------------------
tab16s <- fread(IN_16S)
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
## 7) Full-sample fit (for coefficient tables + nomogram)
## -------------------------
db_m1 <- copy(ps)
db_m1[, z_abund := if (sd(log_Enterococcus_16S) == 0) 0 else as.numeric(scale(log_Enterococcus_16S))]
fit_m1 <- fit_firth_stable(MR ~ z_abund, db_m1)
coef_m1 <- coef_table_from_fit(fit_m1, "M1_Enterococcus_16S")
fwrite(coef_m1, file.path(OUT_DIR, "03_model_coefficients", "M1_Enterococcus_16S_coefficients.csv"))

db_m3 <- copy(ps)
db_m3[, z_abund := if (sd(log_Enterococcus_16S) == 0) 0 else as.numeric(scale(log_Enterococcus_16S))]
db_m3[, z_phage := if (sd(log_Enterococcus_phage) == 0) 0 else as.numeric(scale(log_Enterococcus_phage))]
fit_m3 <- fit_firth_stable(MR ~ z_abund + z_phage, db_m3)
coef_m3 <- coef_table_from_fit(fit_m3, "M3_Combined")
fwrite(coef_m3, file.path(OUT_DIR, "03_model_coefficients", "M3_Combined_coefficients.csv"))

coef_all <- rbindlist(list(coef_m3, coef_m1), fill = TRUE)
fwrite(coef_all, file.path(OUT_DIR, "03_model_coefficients", "M3_M1_coefficients_combined.csv"))

## apparent predictions
ps[, p_app_M1 := clamp_prob(as.numeric(predict(fit_m1, newdata = db_m1, type = "response")))]
ps[, p_app_M3 := clamp_prob(as.numeric(predict(fit_m3, newdata = db_m3, type = "response")))]

## -------------------------
## 8) LOOCV prediction
## -------------------------
message("LOOCV: M1_Enterococcus_16S")
ps[, p_loocv_M1 := make_loocv_pred("M1_Enterococcus_16S", ps)]
message("LOOCV: M3_Combined")
ps[, p_loocv_M3 := make_loocv_pred("M3_Combined", ps)]

pred_out <- ps[, .(
  Sample, Group, MR,
  Enterococcus_16S_counts, log_Enterococcus_16S,
  Enterococcus_phage_burden_TPM, log_Enterococcus_phage,
  p_app_M1, p_loocv_M1,
  p_app_M3, p_loocv_M3
)]
fwrite(pred_out, file.path(OUT_DIR, "04_loocv_predictions", "M3_primary_M1_support_per_sample_predictions.csv"))

## -------------------------
## 9) Metrics / ROC points / best thresholds
## -------------------------
roc_m1 <- roc_points(ps$MR, ps$p_loocv_M1)
roc_m3 <- roc_points(ps$MR, ps$p_loocv_M3)
auc_m1 <- auc_trapz(roc_m1)
auc_m3 <- auc_trapz(roc_m3)

bt_m1 <- best_threshold_youden(ps$MR, ps$p_loocv_M1)
bt_m3 <- best_threshold_youden(ps$MR, ps$p_loocv_M3)

metric_tab <- rbindlist(list(
  cbind(data.table(model = "M3_Combined", pretty_label = PRETTY_MAP[["M3_Combined"]], AUC_LOOCV = auc_m3), bt_m3),
  cbind(data.table(model = "M1_Enterococcus_16S", pretty_label = PRETTY_MAP[["M1_Enterococcus_16S"]], AUC_LOOCV = auc_m1), bt_m1)
), fill = TRUE)
metric_tab[, preferred_role := c("Primary", "Support")]

fwrite(metric_tab, file.path(OUT_DIR, "05_metrics", "M3_M1_LOOCV_AUC_threshold_summary.csv"))

rp_m1 <- copy(roc_m1); rp_m1[, `:=`(model = "M1_Enterococcus_16S", pretty_label = PRETTY_MAP[["M1_Enterococcus_16S"]])]
rp_m3 <- copy(roc_m3); rp_m3[, `:=`(model = "M3_Combined", pretty_label = PRETTY_MAP[["M3_Combined"]])]
fwrite(rp_m1, file.path(OUT_DIR, "06_roc_points", "roc_points_M1_Enterococcus_16S.csv"))
fwrite(rp_m3, file.path(OUT_DIR, "06_roc_points", "roc_points_M3_Combined.csv"))

## -------------------------
## 10) Calibration tables / stats
## -------------------------
cal_m1 <- make_calibration_table(ps$MR, ps$p_loocv_M1, "M1_Enterococcus_16S")
cal_m3 <- make_calibration_table(ps$MR, ps$p_loocv_M3, "M3_Combined")

cal_stats_m1 <- calibration_stats(ps$MR, ps$p_loocv_M1, "M1_Enterococcus_16S")
cal_stats_m3 <- calibration_stats(ps$MR, ps$p_loocv_M3, "M3_Combined")
cal_stats_all <- rbindlist(list(cal_stats_m3, cal_stats_m1), fill = TRUE)

fwrite(cal_m1, file.path(OUT_DIR, "07_calibration_tables", "calibration_bin_table_M1.csv"))
fwrite(cal_m3, file.path(OUT_DIR, "07_calibration_tables", "calibration_bin_table_M3.csv"))
fwrite(cal_stats_all, file.path(OUT_DIR, "07_calibration_tables", "calibration_stats_M3_M1.csv"))

## -------------------------
## 11) DCA tables
## -------------------------
dca_m1 <- make_dca_table(ps$MR, ps$p_loocv_M1, "M1_Enterococcus_16S")
dca_m3 <- make_dca_table(ps$MR, ps$p_loocv_M3, "M3_Combined")

fwrite(dca_m1, file.path(OUT_DIR, "08_dca_tables", "DCA_threshold_table_M1.csv"))
fwrite(dca_m3, file.path(OUT_DIR, "08_dca_tables", "DCA_threshold_table_M3.csv"))

## -------------------------
## 12) Nomogram point tables
## -------------------------
nom_m1 <- make_nomogram_axes(fit_m1, db_m1, "M1_Enterococcus_16S")
nom_m3 <- make_nomogram_axes(fit_m3, db_m3, "M3_Combined")

fwrite(nom_m1$ax_ab, file.path(OUT_DIR, "09_nomogram_points", "Nomogram_axis_M1_z_abund.csv"))
fwrite(nom_m1$ax_total, file.path(OUT_DIR, "09_nomogram_points", "Nomogram_totalpoints_to_prob_M1.csv"))

fwrite(nom_m3$ax_ab, file.path(OUT_DIR, "09_nomogram_points", "Nomogram_axis_M3_z_abund.csv"))
fwrite(nom_m3$ax_ph, file.path(OUT_DIR, "09_nomogram_points", "Nomogram_axis_M3_z_phage.csv"))
fwrite(nom_m3$ax_total, file.path(OUT_DIR, "09_nomogram_points", "Nomogram_totalpoints_to_prob_M3.csv"))

## per-sample nomogram tables
nd_m1 <- copy(db_m1)
nd_m1[, points_abund := nom_m1$beta_ab * z_abund * nom_m1$scale_points]
nd_m1[, total_points := points_abund]
nd_m1[, p_from_points := plogis(nom_m1$intercept + total_points / nom_m1$scale_points)]
fwrite(nd_m1[, .(Sample = ps$Sample, Group = ps$Group, MR = ps$MR, log_Enterococcus_16S = ps$log_Enterococcus_16S,
                 z_abund, points_abund, total_points, p_from_points, p_app_M1 = ps$p_app_M1, p_loocv_M1 = ps$p_loocv_M1)],
       file.path(OUT_DIR, "09_nomogram_points", "Nomogram_points_per_sample_M1.csv"))

nd_m3 <- copy(db_m3)
nd_m3[, points_abund := nom_m3$beta_ab * z_abund * nom_m3$scale_points]
nd_m3[, points_phage := nom_m3$beta_ph * z_phage * nom_m3$scale_points]
nd_m3[, total_points := points_abund + points_phage]
nd_m3[, p_from_points := plogis(nom_m3$intercept + total_points / nom_m3$scale_points)]
fwrite(nd_m3[, .(Sample = ps$Sample, Group = ps$Group, MR = ps$MR,
                 log_Enterococcus_16S = ps$log_Enterococcus_16S,
                 log_Enterococcus_phage = ps$log_Enterococcus_phage,
                 z_abund, z_phage,
                 points_abund, points_phage, total_points,
                 p_from_points, p_app_M3 = ps$p_app_M3, p_loocv_M3 = ps$p_loocv_M3)],
       file.path(OUT_DIR, "09_nomogram_points", "Nomogram_points_per_sample_M3.csv"))

## -------------------------
## 13) Plots: M3 primary
## -------------------------
plot_roc_single(
  roc_dt = roc_m3,
  auc_val = auc_m3,
  title_txt = "M3 ROC (LOOCV): abundance + phage burden",
  color_hex = MODEL_COLORS[["M3_Combined"]],
  out_pdf = file.path(OUT_DIR, "10_plots_M3_primary", "ROC_M3_Combined.pdf"),
  out_png = file.path(OUT_DIR, "10_plots_M3_primary", "ROC_M3_Combined.png")
)

plot_calibration_single(
  cal_tab = cal_m3,
  stats_dt = cal_stats_m3,
  title_txt = "Calibration (LOOCV): M3 primary model",
  color_hex = MODEL_COLORS[["M3_Combined"]],
  out_pdf = file.path(OUT_DIR, "10_plots_M3_primary", "Calibration_M3_Combined.pdf"),
  out_png = file.path(OUT_DIR, "10_plots_M3_primary", "Calibration_M3_Combined.png")
)

plot_dca_single(
  dca_tab = dca_m3,
  title_txt = "Decision Curve Analysis (LOOCV): M3 primary model",
  model_color = MODEL_COLORS[["M3_Combined"]],
  out_pdf = file.path(OUT_DIR, "10_plots_M3_primary", "DCA_M3_Combined.pdf"),
  out_png = file.path(OUT_DIR, "10_plots_M3_primary", "DCA_M3_Combined.png")
)

plot_nomogram_axis(
  ax_dt = nom_m3$ax_ab,
  title_txt = "Nomogram component (M3): z_Enterococcus abundance",
  xlab_txt = "z_Enterococcus abundance",
  out_pdf = file.path(OUT_DIR, "10_plots_M3_primary", "Nomogram_component_M3_z_abund.pdf"),
  out_png = file.path(OUT_DIR, "10_plots_M3_primary", "Nomogram_component_M3_z_abund.png")
)
plot_nomogram_axis(
  ax_dt = nom_m3$ax_ph,
  title_txt = "Nomogram component (M3): z_Enterococcus phage burden",
  xlab_txt = "z_Enterococcus phage burden",
  out_pdf = file.path(OUT_DIR, "10_plots_M3_primary", "Nomogram_component_M3_z_phage.pdf"),
  out_png = file.path(OUT_DIR, "10_plots_M3_primary", "Nomogram_component_M3_z_phage.png")
)
plot_nomogram_total(
  ax_dt = nom_m3$ax_total,
  title_txt = "Nomogram total points -> predicted probability (M3)",
  out_pdf = file.path(OUT_DIR, "10_plots_M3_primary", "Nomogram_totalpoints_to_prob_M3.pdf"),
  out_png = file.path(OUT_DIR, "10_plots_M3_primary", "Nomogram_totalpoints_to_prob_M3.png")
)

## -------------------------
## 14) Plots: M1 support
## -------------------------
plot_roc_single(
  roc_dt = roc_m1,
  auc_val = auc_m1,
  title_txt = "M1 ROC (LOOCV): Enterococcus genus abundance",
  color_hex = MODEL_COLORS[["M1_Enterococcus_16S"]],
  out_pdf = file.path(OUT_DIR, "11_plots_M1_support", "ROC_M1_Enterococcus_16S.pdf"),
  out_png = file.path(OUT_DIR, "11_plots_M1_support", "ROC_M1_Enterococcus_16S.png")
)

plot_calibration_single(
  cal_tab = cal_m1,
  stats_dt = cal_stats_m1,
  title_txt = "Calibration (LOOCV): M1 support model",
  color_hex = MODEL_COLORS[["M1_Enterococcus_16S"]],
  out_pdf = file.path(OUT_DIR, "11_plots_M1_support", "Calibration_M1_Enterococcus_16S.pdf"),
  out_png = file.path(OUT_DIR, "11_plots_M1_support", "Calibration_M1_Enterococcus_16S.png")
)

plot_dca_single(
  dca_tab = dca_m1,
  title_txt = "Decision Curve Analysis (LOOCV): M1 support model",
  model_color = MODEL_COLORS[["M1_Enterococcus_16S"]],
  out_pdf = file.path(OUT_DIR, "11_plots_M1_support", "DCA_M1_Enterococcus_16S.pdf"),
  out_png = file.path(OUT_DIR, "11_plots_M1_support", "DCA_M1_Enterococcus_16S.png")
)

plot_nomogram_axis(
  ax_dt = nom_m1$ax_ab,
  title_txt = "Nomogram component (M1): z_Enterococcus abundance",
  xlab_txt = "z_Enterococcus abundance",
  out_pdf = file.path(OUT_DIR, "11_plots_M1_support", "Nomogram_component_M1_z_abund.pdf"),
  out_png = file.path(OUT_DIR, "11_plots_M1_support", "Nomogram_component_M1_z_abund.png")
)
plot_nomogram_total(
  ax_dt = nom_m1$ax_total,
  title_txt = "Nomogram total points -> predicted probability (M1)",
  out_pdf = file.path(OUT_DIR, "11_plots_M1_support", "Nomogram_totalpoints_to_prob_M1.pdf"),
  out_png = file.path(OUT_DIR, "11_plots_M1_support", "Nomogram_totalpoints_to_prob_M1.png")
)

## -------------------------
## 15) Combined plots: M3 vs M1
## -------------------------
plot_roc_combined(
  roc_m1 = roc_m1,
  roc_m3 = roc_m3,
  auc_m1 = auc_m1,
  auc_m3 = auc_m3,
  out_pdf = file.path(OUT_DIR, "12_plots_M3_vs_M1", "ROC_M3_vs_M1_combined.pdf"),
  out_png = file.path(OUT_DIR, "12_plots_M3_vs_M1", "ROC_M3_vs_M1_combined.png")
)

plot_calibration_combined(
  cal_m1 = cal_m1,
  cal_m3 = cal_m3,
  out_pdf = file.path(OUT_DIR, "12_plots_M3_vs_M1", "Calibration_M3_vs_M1_combined.pdf"),
  out_png = file.path(OUT_DIR, "12_plots_M3_vs_M1", "Calibration_M3_vs_M1_combined.png")
)

plot_dca_combined(
  dca_m1 = dca_m1,
  dca_m3 = dca_m3,
  out_pdf = file.path(OUT_DIR, "12_plots_M3_vs_M1", "DCA_M3_vs_M1_combined.pdf"),
  out_png = file.path(OUT_DIR, "12_plots_M3_vs_M1", "DCA_M3_vs_M1_combined.png")
)

## -------------------------
## 16) Summary text
## -------------------------
summary_txt <- c(
  paste0("Output directory: ", OUT_DIR),
  "",
  "Follow-up evaluation design:",
  "- Primary model: M3 = Enterococcus genus abundance + phage burden",
  "- Support model: M1 = Enterococcus genus abundance only",
  "- No covariate adjustment",
  "- LOOCV-based predicted probabilities used for ROC / calibration / DCA",
  "- Nomogram based on full-sample Firth logistic fits with standardized predictors",
  "",
  "Key files:",
  "- 02_per_sample_table/per_sample_Genus_unadjusted_16S_plus_phage_table.csv",
  "- 04_loocv_predictions/M3_primary_M1_support_per_sample_predictions.csv",
  "- 05_metrics/M3_M1_LOOCV_AUC_threshold_summary.csv",
  "- 07_calibration_tables/calibration_stats_M3_M1.csv",
  "- 10_plots_M3_primary/*",
  "- 11_plots_M1_support/*",
  "- 12_plots_M3_vs_M1/*",
  "",
  paste0("M3 AUC_LOOCV: ", sprintf("%.4f", auc_m3)),
  paste0("M1 AUC_LOOCV: ", sprintf("%.4f", auc_m1))
)
writeLines(summary_txt, con = file.path(OUT_DIR, "README_summary.txt"))

message("\n✅ 完成：M3 主分析 + M1 补充性对照分析已输出。")
message("优先查看：")
message("  1) 05_metrics/M3_M1_LOOCV_AUC_threshold_summary.csv")
message("  2) 10_plots_M3_primary/ROC_M3_Combined.png")
message("  3) 10_plots_M3_primary/Calibration_M3_Combined.png")
message("  4) 10_plots_M3_primary/DCA_M3_Combined.png")
message("  5) 10_plots_M3_primary/Nomogram_totalpoints_to_prob_M3.png")
message("  6) 11_plots_M1_support/ROC_M1_Enterococcus_16S.png")
message("  7) 12_plots_M3_vs_M1/ROC_M3_vs_M1_combined.png")
message("  8) 12_plots_M3_vs_M1/Calibration_M3_vs_M1_combined.png")
message("  9) 12_plots_M3_vs_M1/DCA_M3_vs_M1_combined.png")
