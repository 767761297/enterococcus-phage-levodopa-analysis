#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

## ================================================================
## 16S–Virome correlation network-heatmap style plot
## Transparent background + linewidth-coded significance + italic phage labels
## FIXED2: add per-cell heatmap borders like reference figure,
##         remove outer frame, restore linewidth legend, keep one rho legend
## ================================================================

## ---------- user options ----------
INPUT_DIR <- NA_character_

P_CUTOFF              <- 0.05
EDGE_CURVATURE        <- 0.15
EDGE_ALPHA            <- 0.72
LEFT_HEATMAP_SIZE_MIN <- 2.2
LEFT_HEATMAP_SIZE_MAX <- 8.6

## significance -> line width legend
SIG_LEVELS <- c(
  "p < 0.05 & q > 0.1",
  "p < 0.05 & q < 0.1",
  "p < 0.05 & q < 0.05"
)
SIG_WIDTHS <- c(
  "p < 0.05 & q > 0.1" = 0.5,
  "p < 0.05 & q < 0.1" = 1.5,
  "p < 0.05 & q < 0.05" = 2.00
)
SIG_LABELS <- c(
  "p < 0.05 & q > 0.1" = "*\np < 0.05 & q > 0.1",
  "p < 0.05 & q < 0.1" = "**\np < 0.05 & q < 0.1",
  "p < 0.05 & q < 0.05" = "***\np < 0.05 & q < 0.05"
)
## layout
HEATMAP_LEFT_PAD      <- 1.9
HEATMAP_BOTTOM_PAD    <- 1.15
ANCHOR_X              <- 15.0
NODE_X                <- 19.9
NODE_LABEL_OFFSET     <- 0.72
RIGHT_X_EXTRA         <- 10.0
TOP_BOTTOM_PAD        <- 0.8

## export sizes
PDF_W                 <- 18.6
PDF_H                 <- 10.6
PNG_RES               <- 300

find_input_dir <- function() {
  candidates <- c(
    "results_16S_virome_bubbleplot",
    "results_16S_virome_selected_genus_correlation_unadjusted_revised_stars",
    "results_16S_virome_selected_genus_correlation_unadjusted",
    list.dirs(".", recursive = FALSE, full.names = FALSE)
  )
  candidates <- unique(candidates[!is.na(candidates) & candidates != ""])

  needed_16  <- file.path("02_per_sample_matrices", "selected_16S_genus_CLR_samples_x_features.tsv")
  needed_vir <- file.path("02_per_sample_matrices", "selected_virome_host_genus_CLR_samples_x_features.tsv")
  long_opts  <- c(
    file.path("03_stats", "pairwise_spearman_CLR_16S_vs_virome_host_genus_long_revised_stars.tsv"),
    file.path("03_stats", "pairwise_spearman_CLR_16S_vs_virome_host_genus_long.tsv")
  )

  ok <- vapply(candidates, function(d) {
    file.exists(file.path(d, needed_16)) &&
      file.exists(file.path(d, needed_vir)) &&
      any(file.exists(file.path(d, long_opts)))
  }, logical(1))

  hits <- candidates[ok]
  if (length(hits) == 0) stop("Could not auto-detect previous output directory. Set INPUT_DIR manually.")
  hits[1]
}

pick_existing <- function(paths) {
  hit <- paths[file.exists(paths)][1]
  if (is.na(hit)) stop("None of these files exists: ", paste(paths, collapse = " | "))
  hit
}

get_sample_col <- function(df, file_label) {
  cand <- c("sample_id", "Sample", "sample", "SampleID")
  hit <- cand[cand %in% names(df)][1]
  if (is.na(hit)) stop(file_label, " must contain one sample column among: ", paste(cand, collapse = ", "))
  hit
}

rescale_range <- function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {
  if (any(!is.finite(from)) || diff(from) == 0) return(rep(mean(to), length(x)))
  (x - from[1]) / diff(from) * diff(to) + to[1]
}

short_virome_label <- function(x) sub("\\s*\\|.*$", "", x)

make_phage_parse_label <- function(x) {
  x <- short_virome_label(x)
  x <- gsub("'", "\\\\'", x, fixed = TRUE)
  sprintf("italic('%s')*'_phage'", x)
}

classify_sig <- function(p, q) {
  out <- rep(NA_character_, length(p))
  out[!is.na(p) & !is.na(q) & p < 0.05 & q > 0.10] <- SIG_LEVELS[1]
  out[!is.na(p) & !is.na(q) & p < 0.05 & q < 0.10] <- SIG_LEVELS[2]
  out[!is.na(p) & !is.na(q) & p < 0.05 & q < 0.05] <- SIG_LEVELS[3]
  factor(out, levels = SIG_LEVELS, ordered = TRUE)
}

calc_pairwise_spearman <- function(mat) {
  feats <- colnames(mat)
  cor_mat <- suppressWarnings(cor(mat, method = "spearman", use = "pairwise.complete.obs"))
  p_mat   <- matrix(NA_real_, nrow = ncol(mat), ncol = ncol(mat), dimnames = list(feats, feats))
  q_mat   <- matrix(NA_real_, nrow = ncol(mat), ncol = ncol(mat), dimnames = list(feats, feats))
  n_mat   <- matrix(NA_integer_, nrow = ncol(mat), ncol = ncol(mat), dimnames = list(feats, feats))

  for (i in seq_along(feats)) {
    for (j in seq_along(feats)) {
      x <- mat[, i]
      y <- mat[, j]
      ok <- complete.cases(x, y)
      n_mat[i, j] <- sum(ok)
      if (sum(ok) >= 3) {
        if (i == j) {
          p_mat[i, j] <- 0
        } else {
          ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
          p_mat[i, j] <- ct$p.value
        }
      }
    }
  }

  lt_vals <- p_mat[lower.tri(p_mat, diag = FALSE)]
  q_vals <- rep(NA_real_, length(lt_vals))
  keep <- !is.na(lt_vals)
  q_vals[keep] <- p.adjust(lt_vals[keep], method = "BH")
  q_mat[lower.tri(q_mat, diag = FALSE)] <- q_vals
  q_mat <- t(q_mat)
  q_mat[lower.tri(q_mat, diag = FALSE)] <- q_vals
  diag(q_mat) <- 0

  long <- expand.grid(feature_y = feats, feature_x = feats, stringsAsFactors = FALSE) |>
    as_tibble() |>
    rowwise() |>
    mutate(
      rho = cor_mat[feature_y, feature_x],
      p_value = p_mat[feature_y, feature_x],
      q_value = q_mat[feature_y, feature_x],
      n_samples = n_mat[feature_y, feature_x]
    ) |>
    ungroup()

  list(cor_mat = cor_mat, p_mat = p_mat, q_mat = q_mat, n_mat = n_mat, long = long)
}

## ---------- locate previous output ----------
if (is.na(INPUT_DIR) || INPUT_DIR == "") INPUT_DIR <- find_input_dir()

f16  <- file.path(INPUT_DIR, "02_per_sample_matrices", "selected_16S_genus_CLR_samples_x_features.tsv")
fvir <- file.path(INPUT_DIR, "02_per_sample_matrices", "selected_virome_host_genus_CLR_samples_x_features.tsv")
fcor <- pick_existing(c(
  file.path(INPUT_DIR, "03_stats", "pairwise_spearman_CLR_16S_vs_virome_host_genus_long_revised_stars.tsv"),
  file.path(INPUT_DIR, "03_stats", "pairwise_spearman_CLR_16S_vs_virome_host_genus_long.tsv")
))
stopifnot(file.exists(f16), file.exists(fvir), file.exists(fcor))

OUTDIR <- file.path(INPUT_DIR, "06_network_heatmap")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "01_stats"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "02_plot_data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "03_figures_pdf"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "04_figures_png"), recursive = TRUE, showWarnings = FALSE)

message("Using previous-output directory: ", normalizePath(INPUT_DIR))
d16  <- fread(f16)
dvir <- fread(fvir)
dcor <- fread(fcor)

sample_col_16  <- get_sample_col(d16,  "16S CLR matrix")
sample_col_vir <- get_sample_col(dvir, "Virome CLR matrix")

samples_common <- intersect(d16[[sample_col_16]], dvir[[sample_col_vir]])
if (length(samples_common) < 3) stop("Too few overlapping samples between 16S and virome matrices.")

common_order <- d16[[sample_col_16]][d16[[sample_col_16]] %in% samples_common]
d16  <- d16  |> filter(.data[[sample_col_16]] %in% common_order) |> arrange(match(.data[[sample_col_16]], common_order))
dvir <- dvir |> filter(.data[[sample_col_vir]] %in% common_order) |> arrange(match(.data[[sample_col_vir]], common_order))

feat16  <- setdiff(names(d16), sample_col_16)
featvir <- setdiff(names(dvir), sample_col_vir)

mat16 <- as.matrix(d16[, ..feat16])
rownames(mat16) <- d16[[sample_col_16]]
storage.mode(mat16) <- "numeric"

## ---------- 16S internal correlation ----------
res16 <- calc_pairwise_spearman(mat16)
write.table(res16$cor_mat, file.path(OUTDIR, "01_stats", "selected_16S_internal_spearman_rho_matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(res16$p_mat,   file.path(OUTDIR, "01_stats", "selected_16S_internal_spearman_p_matrix.tsv"),   sep = "\t", quote = FALSE, col.names = NA)
write.table(res16$q_mat,   file.path(OUTDIR, "01_stats", "selected_16S_internal_spearman_q_matrix.tsv"),   sep = "\t", quote = FALSE, col.names = NA)
write.table(res16$long,    file.path(OUTDIR, "01_stats", "selected_16S_internal_spearman_long.tsv"),       sep = "\t", quote = FALSE, row.names = FALSE)

n16 <- length(feat16)
feat16_y <- rev(feat16)
map_y <- setNames(seq_len(n16), feat16_y)
map_x <- setNames(seq_len(n16), feat16)

lower_df <- expand.grid(feature_y = feat16, feature_x = feat16, stringsAsFactors = FALSE) |>
  as_tibble() |>
  mutate(i = match(feature_y, feat16), j = match(feature_x, feat16)) |>
  filter(i >= j) |>
  mutate(
    rho = mapply(function(a, b) res16$cor_mat[a, b], feature_y, feature_x),
    p_value = mapply(function(a, b) res16$p_mat[a, b], feature_y, feature_x),
    q_value = mapply(function(a, b) res16$q_mat[a, b], feature_y, feature_x),
    x = unname(map_x[feature_x]),
    y = unname(map_y[feature_y]),
    abs_rho = abs(rho),
    star = case_when(
      i == j ~ "",
      p_value < 0.05 & q_value < 0.05 ~ "***",
      p_value < 0.05 & q_value < 0.10 ~ "**",
      p_value < 0.05 & q_value > 0.10 ~ "*",
      TRUE ~ ""
    )
  )
write.table(lower_df, file.path(OUTDIR, "02_plot_data", "left_heatmap_lower_triangle_plotdata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
lower_star_df <- lower_df |>
  filter(i > j, star != "")
write.table(lower_star_df, file.path(OUTDIR, "02_plot_data", "left_heatmap_lower_triangle_stars.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

lower_grid_df <- expand.grid(feature_y = feat16, feature_x = feat16, stringsAsFactors = FALSE) |>
  as_tibble() |>
  mutate(i = match(feature_y, feat16), j = match(feature_x, feat16)) |>
  filter(i >= j) |>
  mutate(
    x = unname(map_x[feature_x]),
    y = unname(map_y[feature_y])
  )
write.table(lower_grid_df, file.path(OUTDIR, "02_plot_data", "left_heatmap_lower_triangle_gridcells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---------- cross-omics edges ----------
need_base <- c("feature_16S_genus", "rho", "p_value", "q_value")
miss <- setdiff(need_base, names(dcor))
if (length(miss) > 0) stop("Cross-omics long table missing columns: ", paste(miss, collapse = ", "))

if ("feature_virome_display" %in% names(dcor)) {
  dcor$feature_virome_display <- as.character(dcor$feature_virome_display)
} else if ("feature_virome_host_genus" %in% names(dcor)) {
  dcor$feature_virome_display <- as.character(dcor$feature_virome_host_genus)
} else {
  stop("Cross-omics long table must contain 'feature_virome_display' or 'feature_virome_host_genus'.")
}

virome_display_order <- featvir
nvir <- length(virome_display_order)
virome_short_order <- short_virome_label(virome_display_order)
virome_parse_order <- vapply(virome_display_order, make_phage_parse_label, character(1))

anchor_y_df <- tibble(
  feature_16S_genus = feat16,
  x_start = ANCHOR_X,
  y_start = unname(map_y[feat16])
)
node_y_df <- tibble(
  feature_virome_display = virome_display_order,
  feature_virome_short   = virome_short_order,
  feature_virome_parse   = virome_parse_order,
  x_end = NODE_X,
  y_end = seq(from = n16, to = 1, length.out = nvir)
)

edges_all <- dcor |>
  mutate(
    feature_16S_genus = as.character(feature_16S_genus),
    feature_virome_display = as.character(feature_virome_display),
    feature_virome_short = short_virome_label(feature_virome_display)
  ) |>
  left_join(anchor_y_df, by = "feature_16S_genus") |>
  left_join(node_y_df, by = c("feature_virome_display", "feature_virome_short")) |>
  filter(!is.na(x_start), !is.na(x_end)) |>
  mutate(
    abs_rho = abs(rho),
    keep_edge = p_value < P_CUTOFF,
    sig_class = classify_sig(p_value, q_value)
  )

edges_plot <- edges_all |>
  filter(keep_edge, !is.na(sig_class))

write.table(edges_all,  file.path(OUTDIR, "02_plot_data", "cross_omics_edges_all_with_coordinates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(edges_plot, file.path(OUTDIR, "02_plot_data", "cross_omics_edges_shown_p_lt_0.05.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

bottom_lab_df <- tibble(feature = feat16, x = seq_len(n16), y = 0.12)
left_lab_df   <- tibble(feature = feat16, x = 0.10, y = unname(map_y[feat16]))

write.table(node_y_df, file.path(OUTDIR, "02_plot_data", "right_virome_node_labels_italicphage.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

subtitle_txt <- paste0(
  "Left: internal Spearman correlation among selected 16S genera; right edges: cross-omics Spearman correlations (p < 0.05 shown)\n",
  "Edge linewidth encodes significance: * p < 0.05 & q > 0.1 / ** p < 0.05 & q < 0.1 / *** p < 0.05 & q < 0.05 (q = FDR)"
)

p <- ggplot() +
  geom_tile(
    data = lower_df,
    aes(x = x, y = y, fill = rho),
    width = 0.82, height = 0.82,
    color = "grey60", linewidth = 0.42
  ) +
  geom_text(
    data = lower_star_df,
    aes(x = x, y = y, label = star),
    size = 4,
    fontface = "bold",
    color = "black"
  ) +
  geom_curve(
    data = edges_plot,
    aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = rho, linewidth = sig_class),
    curvature = EDGE_CURVATURE,
    alpha = EDGE_ALPHA,
    ncp = 20,
    lineend = "round",
    show.legend = TRUE
  ) +
  geom_point(
    data = anchor_y_df,
    aes(x = x_start, y = y_start),
    shape = 21, size = 2.5, stroke = 0.55,
    fill = "#f06292", color = "#3f51b5"
  ) +
  geom_point(
    data = node_y_df,
    aes(x = x_end, y = y_end),
    shape = 21, size = 3.0, stroke = 0.7,
    fill = "#f06292", color = "#3f51b5"
  ) +
  geom_text(
    data = node_y_df,
    aes(x = x_end + NODE_LABEL_OFFSET, y = y_end, label = feature_virome_parse),
    hjust = 0, size = 4.1, color = "black", parse = TRUE
  ) +
  geom_text(
    data = bottom_lab_df,
    aes(x = x, y = y, label = feature),
    angle = 90, vjust = 0.5, hjust = 1, size = 3.6, color = "black"
  ) +
  geom_text(
    data = left_lab_df,
    aes(x = x, y = y, label = feature),
    hjust = 1, size = 3.6, color = "black"
  ) +
  scale_fill_gradient2(
    low = "#3B6BD6", mid = "white", high = "#D81B36", midpoint = 0,
    limits = c(-1, 1), guide = "none"
  ) +
  scale_color_gradient2(
    low = "#3B6BD6", mid = "white", high = "#D81B36", midpoint = 0,
    limits = c(-1, 1), name = "Spearman's r"
  ) +
  scale_linewidth_manual(
    values = SIG_WIDTHS,
    limits = SIG_LEVELS,
    breaks = SIG_LEVELS,
    labels = SIG_LABELS,
    drop = FALSE,
    name = "Edge significance"
  ) +
  coord_cartesian(
    xlim = c(-HEATMAP_LEFT_PAD, NODE_X + RIGHT_X_EXTRA),
    ylim = c(1 - TOP_BOTTOM_PAD, n16 + TOP_BOTTOM_PAD),
    clip = "off"
  ) +
  labs(
    title = "Selected 16S genus vs virome host-genus correlation network heatmap",
    subtitle = subtitle_txt
  ) +
  guides(
    color = guide_colorbar(order = 1),
    linewidth = guide_legend(
      order = 2,
      override.aes = list(
        color = rep("grey35", length(SIG_LEVELS)),
        alpha = rep(1, length(SIG_LEVELS)),
        linewidth = unname(SIG_WIDTHS),
        linetype = rep(1, length(SIG_LEVELS))
      )
    )
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.background   = element_rect(fill = "transparent", color = NA),
    panel.background  = element_rect(fill = "transparent", color = NA),
    panel.border      = element_blank(),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = NA, color = NA),
    legend.key        = element_rect(fill = "transparent", color = NA),
    plot.margin       = margin(15, 210, 44, 62),
    legend.position   = c(0.93, 0.48),
    legend.justification = c(0, 0.5),
    legend.box        = "vertical",
    legend.spacing.y  = unit(8, "pt"),
    legend.key.width  = unit(18, "pt"),
    legend.key.height = unit(15, "pt"),
    legend.title      = element_text(size = 10.5),
    legend.text       = element_text(size = 9.2, lineheight = 0.95),
    plot.title        = element_text(face = "bold", size = 14, color = "black"),
    plot.subtitle     = element_text(size = 10.1, color = "black")
  )

ggsave(file.path(OUTDIR, "03_figures_pdf", "network_heatmap.pdf"), p, width = PDF_W, height = PDF_H, bg = "transparent")
ggsave(file.path(OUTDIR, "04_figures_png", "network_heatmap.png"), p, width = PDF_W, height = PDF_H, dpi = PNG_RES, bg = "transparent")

cat("========================================\n")
cat("Done.\n")
cat("Input directory used:\n")
cat(normalizePath(INPUT_DIR), "\n")
cat("Output directory:\n")
cat(normalizePath(OUTDIR), "\n")
cat("========================================\n")
