## =========================================
## Enterococcus host-species phage alpha diversity: single-metric notch boxplots (base R)
## Final template aligned to user's approved 16S alpha style
## =========================================

options(stringsAsFactors = FALSE)

suppressWarnings({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Please install it first: install.packages('data.table')")
  }
})

library(data.table)

## ---------- user settings ----------
input_file <- "alpha_per_sample_EnterococcusHostSpeciesPhage.tsv"
root_outdir <- "results_alpha_diversity"
outdir_name <- "02_Enterococcus_host_species_phage"

metric_cols <- c(
  "observed_gt0",
  "observed_gt1",
  "shannon",
  "simpson",
  "inv_simpson",
  "exp_shannon",
  "pielou_evenness"
)

metric_labels <- c(
  observed_gt0 = "Observed richness (>0)",
  observed_gt1 = "Observed richness (>1)",
  shannon = "Shannon diversity",
  simpson = "Simpson index",
  inv_simpson = "Inverse Simpson index",
  exp_shannon = "Exponential Shannon",
  pielou_evenness = "Pielou evenness"
)

group_col <- "group"
group_levels <- c("PD_GR", "PD_MR")

## colors
col_border <- c(PD_GR = "#7CAE00", PD_MR = "#F8766D")
col_fill <- c(
  PD_GR = grDevices::adjustcolor("#7CAE00", alpha.f = 0.50),
  PD_MR = grDevices::adjustcolor("#F8766D", alpha.f = 0.50)
)

pdf_width <- 3.5
pdf_height <- 4.5
png_width <- 1200
png_height <- 1500
png_res <- 320

## box style
boxwex_value <- 0.03
staplewex_value <- 0.22
axis_lwd <- 1.8
box_lwd <- 1.8
sig_lwd <- 1.5
range_value <- 1.5
group_at <- c(1.00, 1.05)
xlim_value <- c(0.95, 1.15)
xaxis_left <- xlim_value[1]
xaxis_right <- 1.08

## text sizes
cex_main <- 1.15
cex_lab <- 1.05
cex_axis <- 0.90
cex_sig <- 1.05

## ---------- helpers ----------
get_sig_label <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 1e-4) return("****")
  if (p < 1e-3) return("***")
  if (p < 1e-2) return("**")
  if (p < 5e-2) return("*")
  "ns"
}

safe_ylim <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(0, 1))
  xr <- range(x)
  if (xr[1] == xr[2]) {
    delta <- ifelse(xr[1] == 0, 1, abs(xr[1]) * 0.15)
    return(c(xr[1] - delta, xr[2] + delta))
  }
  pad <- diff(xr) * 0.20
  c(xr[1] - pad * 0.10, xr[2] + pad)
}

add_sig_bracket_base <- function(x1, x2, y, label, frac = 0.05, lwd = 1.5, cex = 1.05) {
  usr <- par("usr")
  yr <- usr[4] - usr[3]
  h <- yr * frac
  segments(x1, y, x1, y - h, lwd = lwd, xpd = NA)
  segments(x1, y, x2, y, lwd = lwd, xpd = NA)
  segments(x2, y, x2, y - h, lwd = lwd, xpd = NA)
  text((x1 + x2) / 2, y + yr * 0.01, labels = label, cex = cex, font = 2, xpd = NA)
}

summary_by_group <- function(dt, metric_name) {
  out <- dt[, .(
    n = sum(is.finite(get(metric_name))),
    mean = mean(get(metric_name), na.rm = TRUE),
    sd = sd(get(metric_name), na.rm = TRUE),
    median = median(get(metric_name), na.rm = TRUE),
    q1 = quantile(get(metric_name), 0.25, na.rm = TRUE, names = FALSE),
    q3 = quantile(get(metric_name), 0.75, na.rm = TRUE, names = FALSE)
  ), by = group]
  out[, metric := metric_name]
  setcolorder(out, c("metric", "group", "n", "mean", "sd", "median", "q1", "q3"))
  out
}

clamp_notch_conf <- function(bp_obj, min_frac_iqr = 0.03) {
  if (is.null(bp_obj$conf) || is.null(dim(bp_obj$conf))) return(bp_obj)

  for (i in seq_len(ncol(bp_obj$stats))) {
    q1 <- bp_obj$stats[2, i]
    q3 <- bp_obj$stats[4, i]
    iqr <- q3 - q1
    if (!is.finite(iqr) || iqr <= 0) {
      eps <- max(abs(q1), abs(q3), 1) * 1e-6
    } else {
      eps <- iqr * min_frac_iqr
    }

    lo <- bp_obj$conf[1, i]
    hi <- bp_obj$conf[2, i]

    lower_bound <- q1 + eps
    upper_bound <- q3 - eps

    if (lower_bound >= upper_bound) {
      mid <- (q1 + q3) / 2
      lower_bound <- mid - eps / 2
      upper_bound <- mid + eps / 2
    }

    lo <- max(lo, lower_bound)
    hi <- min(hi, upper_bound)

    if (lo >= hi) {
      mid <- (q1 + q3) / 2
      lo <- mid - eps / 2
      hi <- mid + eps / 2
    }

    bp_obj$conf[1, i] <- lo
    bp_obj$conf[2, i] <- hi
  }
  bp_obj
}

make_long_dt <- function(dt, metric_cols, metric_labels) {
  id_cols <- setdiff(names(dt), metric_cols)
  parts <- lapply(metric_cols, function(m) {
    tmp <- copy(dt[, ..id_cols])
    tmp[, metric := m]
    tmp[, value := as.numeric(dt[[m]])]
    tmp[, metric_label := metric_labels[[m]]]
    tmp
  })
  rbindlist(parts, use.names = TRUE, fill = TRUE)
}

## ---------- read input ----------
cat("[1/5] Reading input:", input_file, "\n")
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

dt <- fread(input_file)
required_cols <- c(group_col, metric_cols)
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

dt <- dt[get(group_col) %in% group_levels]
dt[, group := factor(get(group_col), levels = group_levels)]
dt[, (metric_cols) := lapply(.SD, as.numeric), .SDcols = metric_cols]

## ---------- output dirs ----------
outdir <- file.path(root_outdir, outdir_name)
pdf_dir <- file.path(outdir, "plots_pdf")
png_dir <- file.path(outdir, "plots_png")
tbl_dir <- file.path(outdir, "tables")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(pdf_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tbl_dir, showWarnings = FALSE, recursive = TRUE)

## ---------- stats ----------
cat("[2/5] Calculating Wilcoxon p-values and summary statistics...\n")
stat_list <- list()
sum_list <- list()
for (metric_name in metric_cols) {
  x1 <- dt[group == group_levels[1], get(metric_name)]
  x2 <- dt[group == group_levels[2], get(metric_name)]
  wt <- tryCatch(
    wilcox.test(x1, x2, exact = FALSE),
    error = function(e) NULL
  )
  pval <- if (is.null(wt)) NA_real_ else wt$p.value
  stat_list[[metric_name]] <- data.table(
    metric = metric_name,
    label = metric_labels[[metric_name]],
    group1 = group_levels[1],
    group2 = group_levels[2],
    n_group1 = sum(is.finite(x1)),
    n_group2 = sum(is.finite(x2)),
    median_group1 = median(x1, na.rm = TRUE),
    median_group2 = median(x2, na.rm = TRUE),
    p_value = pval,
    sig_label = get_sig_label(pval)
  )
  sum_list[[metric_name]] <- summary_by_group(dt, metric_name)
}
stat_dt <- rbindlist(stat_list, use.names = TRUE)
summary_dt <- rbindlist(sum_list, use.names = TRUE)

fwrite(stat_dt, file.path(tbl_dir, "alpha_Enterococcus_host_species_phage_single_metric_wilcox_stats.tsv"), sep = "\t")
fwrite(summary_dt, file.path(tbl_dir, "alpha_Enterococcus_host_species_phage_single_metric_group_summary.tsv"), sep = "\t")

## ---------- long table ----------
cat("[3/5] Exporting long-format table...\n")
long_dt <- make_long_dt(dt, metric_cols, metric_labels)
fwrite(long_dt, file.path(tbl_dir, "alpha_Enterococcus_host_species_phage_per_sample_long.tsv"), sep = "\t")

## ---------- plotting ----------
cat("[4/5] Generating plots...\n")
plot_one_metric <- function(metric_name, device = c("pdf", "png"), outfile) {
  device <- match.arg(device)
  x1 <- dt[group == group_levels[1], get(metric_name)]
  x2 <- dt[group == group_levels[2], get(metric_name)]
  y_all <- c(x1, x2)
  yl <- safe_ylim(y_all)
  yr <- yl[2] - yl[1]
  sig_y <- max(y_all, na.rm = TRUE) + yr * 0.09
  yl[2] <- max(yl[2], sig_y + yr * 0.08)
  sig_lab <- stat_dt[metric == metric_name, sig_label][1]

  raw_bp <- boxplot(
    list(x1, x2),
    plot = FALSE,
    notch = TRUE,
    outline = FALSE,
    range = range_value
  )
  bp <- clamp_notch_conf(raw_bp, min_frac_iqr = 0.03)

  if (device == "pdf") {
    pdf(outfile, width = pdf_width, height = pdf_height, onefile = FALSE, useDingbats = FALSE)
  } else {
    png(outfile, width = png_width, height = png_height, res = png_res, bg = "white")
  }

  op <- par(
    mar = c(4.8, 4.6, 1.8, 0.8),
    mgp = c(2.1, 0.65, 0),
    xaxs = "i",
    yaxs = "i",
    bg = "white",
    family = "sans"
  )
  on.exit({par(op); dev.off()}, add = TRUE)

  bxp(
    bp,
    notch = TRUE,
    outline = FALSE,
    axes = FALSE,
    ann = FALSE,
    xaxt = "n",
    yaxt = "n",
    frame.plot = FALSE,
    at = group_at,
    xlim = xlim_value,
    boxwex = boxwex_value,
    staplewex = staplewex_value,
    whisklty = 1,
    staplelty = 1,
    medlwd = box_lwd,
    whisklwd = box_lwd,
    staplelwd = box_lwd,
    lwd = box_lwd,
    ylim = yl,
    pars = list(
      boxfill = unname(col_fill[group_levels]),
      boxcol = unname(col_border[group_levels]),
      medcol = unname(col_border[group_levels]),
      whiskcol = unname(col_border[group_levels]),
      staplecol = unname(col_border[group_levels]),
      medlwd = box_lwd,
      whisklwd = box_lwd,
      staplelwd = box_lwd
    )
  )

  axis(2, lwd = axis_lwd, lwd.ticks = axis_lwd, cex.axis = cex_axis, las = 1)
  axis(1, at = group_at, labels = group_levels, lwd = axis_lwd, lwd.ticks = axis_lwd, cex.axis = cex_axis)
  mtext(metric_labels[[metric_name]], side = 3, line = 0.5, cex = cex_main, font = 2)
  mtext(metric_labels[[metric_name]], side = 2, line = 2.5, cex = cex_lab, font = 3)

  box(bty = "n")
  segments(
    x0 = xaxis_left, y0 = par("usr")[3],
    x1 = xaxis_right, y1 = par("usr")[3],
    lwd = axis_lwd, xpd = NA
  )
  segments(
    x0 = par("usr")[1], y0 = par("usr")[3],
    x1 = par("usr")[1], y1 = par("usr")[4],
    lwd = axis_lwd, xpd = NA
  )
  add_sig_bracket_base(group_at[1], group_at[2], y = sig_y, label = sig_lab, frac = 0.045, lwd = sig_lwd, cex = cex_sig)
}

for (metric_name in metric_cols) {
  pdf_file <- file.path(pdf_dir, paste0(metric_name, ".pdf"))
  png_file <- file.path(png_dir, paste0(metric_name, ".png"))
  plot_one_metric(metric_name, device = "pdf", outfile = pdf_file)
  plot_one_metric(metric_name, device = "png", outfile = png_file)
}

cat("[5/5] Done.\n")
cat("Output directory:", outdir, "\n")
cat("Main stats table:", file.path(tbl_dir, "alpha_Enterococcus_host_species_phage_single_metric_wilcox_stats.tsv"), "\n")
cat("PDF plots:", pdf_dir, "\n")
cat("PNG plots:", png_dir, "\n")