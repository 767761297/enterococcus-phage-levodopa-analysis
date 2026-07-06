#!/usr/bin/env Rscript
## =========================================================
## Cross-omics alpha-alpha correlation (custom simplified)
## X = Enterococcus-host phage alpha (virome; Host_species level)
## Y = Enterococcus alpha (16S; within-genus)
##
## Custom changes based on user request:
## - only unadjusted analysis
## - heatmap restricted to 3 x 3 alpha-diversity metrics
##   16S   : chao1, shannon, pielou_evenness
##   Virome: observed_gt0, shannon, pielou_evenness
## - scatter plots: ALL 3 x 3 pairwise overall fits only
## - no groupwise scatter plots
## - polished heatmap with soft blue-white-red palette
##
## Outputs:
## - tables: rho / p / q / n matrices + long table
## - figures: 3 x 3 heatmap
## - figures: 9 overall scatter plots
## - per_sample: data used for each scatter plot
## =========================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(janitor)
  library(ggplot2)
})

## ------------------------
## 0) USER SETTINGS
## ------------------------
INPUT_16S_ALPHA    <- "alpha_per_sample_Enterococcus_species.tsv"
INPUT_VIROME_ALPHA <- "alpha_per_sample_EnterococcusHostSpeciesPhage.tsv"
META_FILE          <- "metadata.csv"

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_ROOT <- paste0("results_alpha_alpha_crossomics_unadjusted_3x3_", timestamp)

GROUP_LEVELS <- c("PD_GR", "PD_MR")
COL_MAP <- c(PD_GR = "#7CAE00", PD_MR = "#F8766D")
LAB_MAP <- c(PD_GR = "GR", PD_MR = "MR")

MAKE_PDF <- TRUE
MAKE_PNG <- TRUE

## requested metrics only
METRICS_16S <- c("chao1", "shannon", "pielou_evenness")
METRICS_VIR <- c("observed_gt0", "shannon", "pielou_evenness")

## significance markers on heatmap
## ** : q < 0.05
## *  : q < 0.10
q_to_star <- function(q) {
  ifelse(is.na(q), "",
         ifelse(q < 0.05, "**",
                ifelse(q < 0.10, "*", "")))
}

metric_lab_16S <- c(
  chao1 = "Chao1",
  shannon = "Shannon",
  pielou_evenness = "Pielou evenness"
)
metric_lab_vir <- c(
  observed_gt0 = "Observed (>0)",
  shannon = "Shannon",
  pielou_evenness = "Pielou evenness"
)

## ------------------------
## 1) Helpers
## ------------------------
dir_init <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(path, "00_manifest"), showWarnings = FALSE)
  dir.create(file.path(path, "tables"), showWarnings = FALSE)
  dir.create(file.path(path, "figures_pdf"), showWarnings = FALSE)
  dir.create(file.path(path, "figures_png"), showWarnings = FALSE)
  dir.create(file.path(path, "per_sample"), showWarnings = FALSE)
}

write_tsv0 <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(x, path)
}

read_tsv_smart <- function(path) {
  tryCatch(
    readr::read_tsv(path, show_col_types = FALSE, locale = locale(encoding = "UTF-8")),
    error = function(e) readr::read_tsv(path, show_col_types = FALSE, locale = locale(encoding = "GB18030"))
  )
}

read_csv_smart <- function(path) {
  tryCatch(
    readr::read_csv(path, show_col_types = FALSE, locale = locale(encoding = "GB18030")),
    error = function(e) readr::read_csv(path, show_col_types = FALSE, locale = locale(encoding = "UTF-8"))
  )
}

normalize_id <- function(x) x |> as.character() |> str_trim() |> str_replace_all("\\s+", "")

find_col <- function(df, patterns) {
  nms <- names(df)
  for (p in patterns) {
    hit <- nms[str_detect(nms, regex(p, ignore_case = TRUE))]
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

spearman_fast <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  n <- length(x)
  if (n < 6) return(c(rho = NA_real_, p = NA_real_, n = n))
  rho <- suppressWarnings(cor(x, y, method = "spearman"))
  if (is.na(rho)) return(c(rho = NA_real_, p = NA_real_, n = n))
  rho2 <- min(0.999999, max(-0.999999, rho))
  tval <- rho2 * sqrt((n - 2) / (1 - rho2^2))
  p <- 2 * pt(-abs(tval), df = n - 2)
  c(rho = rho, p = p, n = n)
}

fmt_p <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.001, formatC(p, format = "e", digits = 2),
                sprintf("%.3f", p)))
}

## ------------------------
## 2) Load metadata
## ------------------------
meta_raw <- read_csv_smart(META_FILE) |> janitor::clean_names()

sample_col <- find_col(meta_raw, c("^sample$"))
group_col  <- find_col(meta_raw, c("^group$"))

stopifnot(!is.na(sample_col), !is.na(group_col))

meta <- meta_raw |>
  mutate(
    sample_id = normalize_id(.data[[sample_col]]),
    group = factor(.data[[group_col]], levels = GROUP_LEVELS)
  ) |>
  filter(group %in% GROUP_LEVELS) |>
  arrange(group, sample_id)

sample_order <- meta$sample_id

## ------------------------
## 3) Load alpha tables
## ------------------------
load_alpha <- function(path, dataset_name) {
  df <- read_tsv_smart(path) |> janitor::clean_names()
  sid <- find_col(df, c("^sample_id$", "^sample$"))
  stopifnot(!is.na(sid))
  df <- df |> mutate(sample_id = normalize_id(.data[[sid]]))

  if (!"group" %in% names(df)) {
    df <- df |> left_join(meta |> select(sample_id, group), by = "sample_id")
  } else {
    df <- df |> mutate(group = factor(as.character(group), levels = GROUP_LEVELS))
  }

  df <- df |>
    filter(sample_id %in% sample_order) |>
    arrange(match(sample_id, sample_order))

  df$dataset <- dataset_name
  df
}

alpha16 <- load_alpha(INPUT_16S_ALPHA, "16S")
alphav  <- load_alpha(INPUT_VIROME_ALPHA, "Virome")

stopifnot(all(METRICS_16S %in% names(alpha16)))
stopifnot(all(METRICS_VIR %in% names(alphav)))

## ------------------------
## 4) Correlation matrices
## ------------------------
run_alpha_alpha_unadjusted <- function(out_dir) {
  dir_init(out_dir)

  M16 <- alpha16 |>
    select(sample_id, all_of(METRICS_16S)) |>
    arrange(match(sample_id, sample_order)) |>
    column_to_rownames("sample_id") |>
    as.matrix()

  MV <- alphav |>
    select(sample_id, all_of(METRICS_VIR)) |>
    arrange(match(sample_id, sample_order)) |>
    column_to_rownames("sample_id") |>
    as.matrix()

  M16[!is.finite(M16)] <- NA_real_
  MV[!is.finite(MV)] <- NA_real_

  rho <- matrix(
    NA_real_,
    nrow = length(METRICS_16S),
    ncol = length(METRICS_VIR),
    dimnames = list(METRICS_16S, METRICS_VIR)
  )
  pmat <- rho
  nmat <- rho

  long_list <- list()
  idx <- 0L

  for (i in seq_along(METRICS_16S)) {
    for (j in seq_along(METRICS_VIR)) {
      st <- spearman_fast(MV[, j], M16[, i])
      idx <- idx + 1L
      long_list[[idx]] <- tibble(
        model = "unadjusted",
        metric_16S = METRICS_16S[i],
        metric_virome = METRICS_VIR[j],
        rho = as.numeric(st["rho"]),
        p_value = as.numeric(st["p"]),
        n_used = as.integer(st["n"])
      )
      rho[METRICS_16S[i], METRICS_VIR[j]] <- as.numeric(st["rho"])
      pmat[METRICS_16S[i], METRICS_VIR[j]] <- as.numeric(st["p"])
      nmat[METRICS_16S[i], METRICS_VIR[j]] <- as.numeric(st["n"])
    }
  }

  long_df <- bind_rows(long_list) |>
    mutate(
      q_value = p.adjust(p_value, method = "BH"),
      star = q_to_star(q_value)
    )

  qmat <- matrix(NA_real_, nrow = nrow(rho), ncol = ncol(rho), dimnames = dimnames(rho))
  for (k in seq_len(nrow(long_df))) {
    qmat[long_df$metric_16S[k], long_df$metric_virome[k]] <- long_df$q_value[k]
  }

  write_tsv0(long_df, file.path(out_dir, "tables", "alpha_alpha_long.tsv"))
  write_tsv0(as.data.frame(rho)  |> rownames_to_column("metric_16S"), file.path(out_dir, "tables", "rho_matrix.tsv"))
  write_tsv0(as.data.frame(pmat) |> rownames_to_column("metric_16S"), file.path(out_dir, "tables", "p_matrix.tsv"))
  write_tsv0(as.data.frame(qmat) |> rownames_to_column("metric_16S"), file.path(out_dir, "tables", "q_matrix.tsv"))
  write_tsv0(as.data.frame(nmat) |> rownames_to_column("metric_16S"), file.path(out_dir, "tables", "n_matrix.tsv"))

  invisible(list(long = long_df, rho = rho, p = pmat, q = qmat, n = nmat))
}

## ------------------------
## 5) Polished heatmap
## ------------------------
plot_heatmap_custom <- function(corr_res, out_dir) {
  df <- corr_res$long |>
    mutate(
      metric_16S_lab = factor(metric_lab_16S[metric_16S], levels = metric_lab_16S[METRICS_16S]),
      metric_virome_lab = factor(metric_lab_vir[metric_virome], levels = metric_lab_vir[METRICS_VIR]),
      rho_label = ifelse(is.na(rho), "NA", sprintf("%.2f", rho)),
      label = ifelse(star == "", rho_label, paste0(rho_label, "\n", star)),
      label_color = ifelse(is.na(rho), "black", ifelse(abs(rho) >= 0.45, "white", "black"))
    )

  p <- ggplot(df, aes(x = metric_virome_lab, y = metric_16S_lab, fill = rho)) +
    geom_tile(width = 0.96, height = 0.96, color = "white", linewidth = 0.9) +
    geom_text(aes(label = label, color = label_color), size = 5.0, fontface = "bold", lineheight = 0.95) +
    scale_color_identity() +
    scale_fill_gradient2(
      low = "#4A77B3",
      mid = "#F7F7F7",
      high = "#C95D63",
      midpoint = 0,
      limits = c(-1, 1),
      breaks = c(-1, -0.5, 0, 0.5, 1),
      name = "Spearman rho"
    ) +
    coord_equal() +
    labs(
      title = "Cross-omics alpha-alpha correlation",
      subtitle = "Virome Enterococcus-host phage vs 16S Enterococcus | unadjusted",
      x = "Virome α-diversity metrics",
      y = "16S α-diversity metrics"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1, color = "black", size = 11),
      axis.text.y = element_text(color = "black", size = 11),
      axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 10)),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11),
      legend.title = element_text(face = "bold"),
      legend.key.height = grid::unit(38, "pt"),
      plot.margin = margin(10, 16, 10, 10)
    )

  if (MAKE_PDF) {
    ggsave(file.path(out_dir, "figures_pdf", "heatmap_rho_unadjusted_3x3.pdf"),
           plot = p, width = 7.2, height = 5.8)
  }
  if (MAKE_PNG) {
    ggsave(file.path(out_dir, "figures_png", "heatmap_rho_unadjusted_3x3.png"),
           plot = p, width = 7.2, height = 5.8, dpi = 320)
  }
}

## ------------------------
## 6) Overall scatter plots only (all 3 x 3 pairs)
## ------------------------
plot_overall_scatter <- function(df, x_col, y_col, x_lab, y_lab, out_pdf, out_png, out_tsv) {
  dat <- df |>
    transmute(
      sample_id,
      group,
      x = .data[[x_col]],
      y = .data[[y_col]]
    ) |>
    filter(is.finite(x), is.finite(y))

  write_tsv0(dat, out_tsv)

  st <- spearman_fast(dat$x, dat$y)
  ann <- paste0(
    "Spearman rho = ", sprintf("%.3f", as.numeric(st["rho"])),
    "\nP = ", fmt_p(as.numeric(st["p"])),
    "\nn = ", as.integer(st["n"])
  )

  p <- ggplot(dat, aes(x = x, y = y)) +
    geom_point(aes(color = group), size = 2.8, alpha = 0.9) +
    geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey75", alpha = 0.28, linewidth = 0.8) +
    scale_color_manual(values = COL_MAP, breaks = names(LAB_MAP), labels = LAB_MAP) +
    labs(
      title = paste0("Overall scatter | ", x_lab, " vs ", y_lab),
      subtitle = "Unadjusted",
      x = x_lab,
      y = y_lab,
      color = "Group"
    ) +
    annotate("text", x = Inf, y = Inf, label = ann, hjust = 1.06, vjust = 1.2, size = 4.1) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10.5),
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.35)
    )

  if (MAKE_PDF) ggsave(out_pdf, p, width = 6.6, height = 5.4)
  if (MAKE_PNG) ggsave(out_png, p, width = 6.6, height = 5.4, dpi = 320)
}

run_all_scatters <- function(out_dir) {
  dir_init(out_dir)

  df <- meta |>
    select(sample_id, group) |>
    left_join(
      alpha16 |>
        select(sample_id, all_of(METRICS_16S)) |>
        rename_with(~ paste0("16S_", .x), all_of(METRICS_16S)),
      by = "sample_id"
    ) |>
    left_join(
      alphav |>
        select(sample_id, all_of(METRICS_VIR)) |>
        rename_with(~ paste0("VIR_", .x), all_of(METRICS_VIR)),
      by = "sample_id"
    ) |>
    arrange(group, sample_id)

  pairs <- expand_grid(
    metric_16S = METRICS_16S,
    metric_vir = METRICS_VIR
  )

  write_tsv0(pairs, file.path(out_dir, "00_manifest", "scatter_pairs_all_3x3.tsv"))

  for (k in seq_len(nrow(pairs))) {
    m16 <- pairs$metric_16S[k]
    mvr <- pairs$metric_vir[k]

    x_col <- paste0("VIR_", mvr)
    y_col <- paste0("16S_", m16)

    x_lab <- paste0("Virome ", metric_lab_vir[[mvr]])
    y_lab <- paste0("16S ", metric_lab_16S[[m16]])

    stem <- paste0("scatter_overall__", x_col, "__", y_col)

    plot_overall_scatter(
      df = df,
      x_col = x_col,
      y_col = y_col,
      x_lab = x_lab,
      y_lab = y_lab,
      out_pdf = file.path(out_dir, "figures_pdf", paste0(stem, ".pdf")),
      out_png = file.path(out_dir, "figures_png", paste0(stem, ".png")),
      out_tsv = file.path(out_dir, "per_sample", paste0(stem, ".tsv"))
    )
  }
}

## ------------------------
## 7) RUN
## ------------------------
dir_init(OUT_ROOT)

write_tsv0(meta, file.path(OUT_ROOT, "00_manifest", "metadata_used.tsv"))
write_tsv0(
  tibble(
    INPUT_16S_ALPHA = INPUT_16S_ALPHA,
    INPUT_VIROME_ALPHA = INPUT_VIROME_ALPHA,
    META_FILE = META_FILE,
    MODEL = "unadjusted"
  ),
  file.path(OUT_ROOT, "00_manifest", "inputs.tsv")
)
write_tsv0(tibble(metric_16S = METRICS_16S, label = unname(metric_lab_16S[METRICS_16S])),
           file.path(OUT_ROOT, "00_manifest", "metrics_16S.tsv"))
write_tsv0(tibble(metric_virome = METRICS_VIR, label = unname(metric_lab_vir[METRICS_VIR])),
           file.path(OUT_ROOT, "00_manifest", "metrics_virome.tsv"))

message("[1/3] Running unadjusted 3x3 alpha-alpha correlations...")
out_a <- file.path(OUT_ROOT, "A_alpha_alpha_unadjusted_3x3")
corr_res <- run_alpha_alpha_unadjusted(out_a)

message("[2/3] Drawing polished heatmap...")
plot_heatmap_custom(corr_res, out_a)

message("[3/3] Drawing 3x3 overall scatter plots...")
out_b <- file.path(OUT_ROOT, "B_scatter_overall_unadjusted_3x3")
run_all_scatters(out_b)

sink(file.path(OUT_ROOT, "00_manifest", "sessionInfo.txt"))
print(sessionInfo())
sink()

message("Done -> ", OUT_ROOT)
