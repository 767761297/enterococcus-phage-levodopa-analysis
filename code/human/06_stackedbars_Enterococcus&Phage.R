#!/usr/bin/env Rscript
# =========================================================
# Stacked bars: 16S Enterococcus species + Enterococcus-host phage (Host_species)
# - Use Source-like input files directly:
#   metadata.csv
#   table_filtered_w_tax.xls   (TSV disguised as .xls; must contain 'taxonomy' + sample columns PDxx)
#   vOTU_abundance_matrix_TPM.tsv
#   cherry_prediction.tsv
#
# Outputs:
#   results_stackedbars_enterococcus_<timestamp>/
#     00_manifest/
#     tables/
#     figures_pdf/
#     figures_png/
# =========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(readr)     # <<< 加这一行
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(cowplot)
  library(jsonlite)
})

# -------------------------
# 0) USER SETTINGS (EDIT IF NEEDED)
# -------------------------
META_FILE   <- "metadata.csv"
TABLE_16S   <- "table_filtered_w_tax.xls"          # TSV disguised as xls
TPM_FILE    <- "vOTU_abundance_matrix_TPM.tsv"
CHERRY_FILE <- "cherry_prediction.tsv"

TOP_N <- 12
STRIP_HEIGHT_RATIO <- 0.045
LEGEND_NCOL <- 2

GROUP_ORDER <- c("PD_GR", "PD_MR")
GROUP_COLORS <- c(PD_GR = "#2ca02c", PD_MR = "#d62728")  # GR green, MR red

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTDIR <- paste0("results_stackedbars_enterococcus_", timestamp)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "00_manifest"), showWarnings = FALSE)
dir.create(file.path(OUTDIR, "tables"), showWarnings = FALSE)
dir.create(file.path(OUTDIR, "figures_pdf"), showWarnings = FALSE)
dir.create(file.path(OUTDIR, "figures_png"), showWarnings = FALSE)

# -------------------------
# 1) Helpers
# -------------------------
is_sample_col <- function(x) str_detect(x, "^PD\\d+$")

parse_tax_rank <- function(tax, rank_char) {
  # tax examples:
  # d_Bacteria;...;g_Enterococcus;s_faecalis
  # d__Bacteria;...;g__Enterococcus;s__faecalis
  if (is.na(tax)) return("")
  parts <- str_split(as.character(tax), ";", simplify = TRUE)
  parts <- str_trim(parts)
  pref1 <- paste0(rank_char, "_")
  pref2 <- paste0(rank_char, "__")
  hit <- parts[str_starts(parts, pref2)]
  if (length(hit) > 0) return(str_remove(hit[1], fixed(pref2)))
  hit <- parts[str_starts(parts, pref1)]
  if (length(hit) > 0) return(str_remove(hit[1], fixed(pref1)))
  ""
}

normalize_species_label <- function(genus, species_val) {
  sp <- str_trim(as.character(species_val))
  if (is.na(sp) || sp == "" || tolower(sp) %in% c("na","nan","unclassified","uncultured","unknown")) {
    return(paste0(genus, "_unclassified"))
  }
  # if it is only epithet (e.g. "faecalis"), add genus
  if (!str_detect(tolower(sp), "^enterococcus")) {
    return(paste0(genus, " ", sp))
  }
  sp
}

# order samples: GR then MR, within group by numeric part of PDxx
order_samples <- function(samples, groups) {
  tibble(sample_id = samples, group = groups) |>
    mutate(
      group = factor(group, levels = GROUP_ORDER),
      sample_num = suppressWarnings(as.integer(str_extract(sample_id, "\\d+")))
    ) |>
    arrange(group, sample_num, sample_id) |>
    pull(sample_id)
}

topn_other <- function(mat, top_n = 12, other_name = "Other") {
  # mat: features x samples, sums to 1 per sample (ideally)
  m <- rowMeans(mat, na.rm = TRUE)
  top_feats <- names(sort(m, decreasing = TRUE))[seq_len(min(top_n, length(m)))]
  out <- mat[top_feats, , drop = FALSE]
  other <- colSums(mat[setdiff(rownames(mat), top_feats), , drop = FALSE], na.rm = TRUE)
  out <- rbind(out, other)
  rownames(out)[nrow(out)] <- other_name
  out
}

mat_to_long <- function(mat, sample_order) {
  as.data.frame(mat) |>
    rownames_to_column("feature") |>
    pivot_longer(-feature, names_to = "sample_id", values_to = "value") |>
    mutate(sample_id = factor(sample_id, levels = sample_order))
}

get_boundary <- function(sample_order, sample_to_group) {
  g <- unname(sample_to_group[sample_order])
  bd <- NA_integer_
  for (i in 2:length(g)) {
    if (!identical(g[i], g[i-1])) { bd <- i; break }
  }
  bd
}

make_strip_plot <- function(sample_order, sample_to_group, boundary_idx) {
  df <- tibble(
    sample_id = factor(sample_order, levels = sample_order),
    group = factor(unname(sample_to_group[sample_order]), levels = GROUP_ORDER)
  )
  p <- ggplot(df, aes(x = sample_id, y = 1, fill = group)) +
    geom_tile(height = 1, width = 1) +
    scale_fill_manual(values = GROUP_COLORS, drop = FALSE) +
    theme_void() +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))
  if (!is.na(boundary_idx)) {
    p <- p + geom_vline(xintercept = boundary_idx - 0.5, linewidth = 0.4, color = "black")
  }
  p
}

make_group_legend <- function() {
  df <- tibble(group = factor(GROUP_ORDER, levels = GROUP_ORDER), x = 1, y = 1)
  p <- ggplot(df, aes(x = x, y = y, fill = group)) +
    geom_tile() +
    scale_fill_manual(values = GROUP_COLORS, drop = FALSE) +
    guides(fill = guide_legend(title = "Group")) +
    theme_void() +
    theme(legend.position = "right")
  cowplot::get_legend(p)
}

make_bar_plot <- function(mat, sample_order, boundary_idx, ylab, title, legend_ncol) {
  df <- mat_to_long(mat, sample_order) |>
    mutate(feature = factor(feature, levels = rownames(mat)[order(rowMeans(mat, na.rm = TRUE), decreasing = TRUE)]))

  p <- ggplot(df, aes(x = sample_id, y = value, fill = feature)) +
    geom_col(width = 0.9, linewidth = 0) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = title,
      x = "Samples (ordered by group; no labels)",
      y = ylab,
      fill = "Features"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(0, 0, 0, 0),
      legend.position = "right"
    ) +
    guides(fill = guide_legend(ncol = legend_ncol))
  if (!is.na(boundary_idx)) {
    p <- p + geom_vline(xintercept = boundary_idx - 0.5, linewidth = 0.4, color = "black")
  }
  p
}

save_plot <- function(plot_obj, stub, n_samples) {
  fig_w <- max(14, n_samples * 0.32)
  fig_h <- 5.2

  ggsave(file.path(OUTDIR, "figures_pdf", paste0(stub, ".pdf")),
         plot_obj, width = fig_w, height = fig_h, device = cairo_pdf)
  ggsave(file.path(OUTDIR, "figures_png", paste0(stub, ".png")),
         plot_obj, width = fig_w, height = fig_h, dpi = 300)
}

plot_stacked_bundle <- function(mat, sample_order, sample_to_group, title, ylab, stub) {
  boundary <- get_boundary(sample_order, sample_to_group)
  strip <- make_strip_plot(sample_order, sample_to_group, boundary)
  bar <- make_bar_plot(mat, sample_order, boundary, ylab, title, LEGEND_NCOL)

  group_leg <- make_group_legend()
  feat_leg <- cowplot::get_legend(bar)

  legends <- cowplot::plot_grid(group_leg, feat_leg, ncol = 1, rel_heights = c(0.22, 0.78))
  bar_noleg <- bar + theme(legend.position = "none")

  main <- cowplot::plot_grid(
    strip, bar_noleg,
    ncol = 1, align = "v",
    rel_heights = c(STRIP_HEIGHT_RATIO, 1 - STRIP_HEIGHT_RATIO)
  )

  final <- cowplot::plot_grid(main, legends, nrow = 1, rel_widths = c(0.73, 0.27))
  save_plot(final, stub, length(sample_order))
}

# -------------------------
# 2) Load metadata -> sample order
# -------------------------
meta <- fread(META_FILE, encoding = "UTF-8")  # if fails on your system, change encoding to "GB18030"
# columns expected: Sample, Group (as in your project)
stopifnot(all(c("Sample", "Group") %in% names(meta)))

meta <- meta |>
  mutate(Sample = as.character(Sample), Group = as.character(Group)) |>
  filter(Group %in% GROUP_ORDER)

sample_order <- order_samples(meta$Sample, meta$Group)
sample_to_group <- setNames(meta$Group, meta$Sample)

write_tsv(meta |> arrange(match(Sample, sample_order)),
          file.path(OUTDIR, "00_manifest", "sample_order_and_group_used.tsv"))

# -------------------------
# 3) 16S: Enterococcus species within-genus RA
# -------------------------
tab16 <- fread(TABLE_16S, sep = "\t", encoding = "UTF-8", data.table = FALSE)
stopifnot("taxonomy" %in% names(tab16))

sample_cols_16 <- names(tab16)[sapply(names(tab16), is_sample_col)]
sample_cols_16 <- intersect(sample_cols_16, sample_order)
stopifnot(length(sample_cols_16) > 0)

tab16$genus <- vapply(tab16$taxonomy, parse_tax_rank, character(1), rank_char = "g")
tab16$species <- vapply(tab16$taxonomy, parse_tax_rank, character(1), rank_char = "s")

ent <- tab16 |>
  filter(tolower(genus) == "enterococcus") |>
  mutate(species_label = vapply(species, function(s) normalize_species_label("Enterococcus", s), character(1))) |>
  select(species_label, all_of(sample_cols_16))

# numeric counts
for (cc in sample_cols_16) ent[[cc]] <- suppressWarnings(as.integer(ent[[cc]]))
ent[is.na(ent)] <- 0L

# aggregate ASVs -> species
sp_counts <- ent |>
  group_by(species_label) |>
  summarise(across(all_of(sample_cols_16), sum), .groups = "drop")

sp_mat <- sp_counts |>
  column_to_rownames("species_label") |>
  as.matrix()
storage.mode(sp_mat) <- "numeric"

# within-genus RA
colsum_ent <- colSums(sp_mat)
ra_16 <- sweep(sp_mat, 2, ifelse(colsum_ent > 0, colsum_ent, NA_real_), "/")
ra_16[is.na(ra_16)] <- 0

write_tsv(as.data.frame(ra_16) |> rownames_to_column("feature"),
          file.path(OUTDIR, "tables", "16S_Enterococcus_species_RA_within_genus.tsv"))

# -------------------------
# 4) Virome: Enterococcus-host phage Host_species RA (suffix _virus)
# -------------------------
tpm <- fread(TPM_FILE, sep = "\t", encoding = "UTF-8", data.table = FALSE)
stopifnot("id" %in% names(tpm))

sample_cols_v <- setdiff(names(tpm), "id")
sample_cols_v <- intersect(sample_cols_v, sample_order)
stopifnot(length(sample_cols_v) > 0)

cherry <- fread(CHERRY_FILE, sep = "\t", encoding = "UTF-8", data.table = FALSE)
stopifnot(all(c("Accession", "Host") %in% names(cherry)))

v <- tpm |>
  left_join(cherry, by = c("id" = "Accession"))

parse_host <- function(host_str) {
  if (is.na(host_str)) return(list(genus = NA_character_, host_species_label = NA_character_))
  s <- str_trim(as.character(host_str))
  m <- str_match(s, "^(species|genus):\\s*(.+)$")
  if (is.na(m[1,1])) return(list(genus = NA_character_, host_species_label = NA_character_))
  rank <- m[1,2]; taxon <- str_trim(m[1,3])
  genus <- str_split(taxon, "\\s+", simplify = TRUE)[1]
  host_species_label <- if (rank == "species") taxon else paste0(genus, " [genus]")
  list(genus = genus, host_species_label = host_species_label)
}

tmp <- lapply(v$Host, parse_host)
v$host_genus <- vapply(tmp, `[[`, character(1), "genus")
v$host_species_label <- vapply(tmp, `[[`, character(1), "host_species_label")

v_ent <- v |>
  filter(host_genus == "Enterococcus") |>
  select(host_species_label, all_of(sample_cols_v))

# numeric TPM
for (cc in sample_cols_v) v_ent[[cc]] <- suppressWarnings(as.numeric(v_ent[[cc]]))
v_ent[is.na(v_ent)] <- 0

hs_tpm <- v_ent |>
  group_by(host_species_label) |>
  summarise(across(all_of(sample_cols_v), sum), .groups = "drop") |>
  mutate(feature = paste0(host_species_label, "_virus")) |>
  select(feature, all_of(sample_cols_v))

hs_mat <- hs_tpm |>
  column_to_rownames("feature") |>
  as.matrix()
storage.mode(hs_mat) <- "numeric"

colsum_v <- colSums(hs_mat)
ra_v <- sweep(hs_mat, 2, ifelse(colsum_v > 0, colsum_v, NA_real_), "/")
ra_v[is.na(ra_v)] <- 0

write_tsv(as.data.frame(ra_v) |> rownames_to_column("feature"),
          file.path(OUTDIR, "tables", "Virome_EnterococcusHostSpeciesPhage_RA_within_subset.tsv"))

# -------------------------
# 5) Plot: all features + TopN+Other (both, same params)
# -------------------------
# All features
plot_stacked_bundle(
  ra_16[, sample_order, drop=FALSE],
  sample_order, sample_to_group,
  title = "16S Enterococcus species composition (within-genus RA; all features)",
  ylab  = "Relative abundance within Enterococcus genus",
  stub  = "stackedbar_16S_Enterococcus_species_all"
)

plot_stacked_bundle(
  ra_v[, sample_order, drop=FALSE],
  sample_order, sample_to_group,
  title = "Enterococcus-host predicted phage composition (Host_species level; all features)",
  ylab  = "Relative abundance within Enterococcus-host phage subset",
  stub  = "stackedbar_Virome_EnterococcusHostSpeciesPhage_all"
)

# TopN + Other
ra_16_top <- topn_other(ra_16[, sample_order, drop=FALSE], top_n = TOP_N, other_name = "Other")
ra_v_top  <- topn_other(ra_v[, sample_order, drop=FALSE],  top_n = TOP_N, other_name = "Other")

write_tsv(tibble(feature = rownames(ra_16_top)),
          file.path(OUTDIR, "tables", "16S_topN_other_features.tsv"))
write_tsv(tibble(feature = rownames(ra_v_top)),
          file.path(OUTDIR, "tables", "Virome_topN_other_features.tsv"))

plot_stacked_bundle(
  ra_16_top,
  sample_order, sample_to_group,
  title = paste0("16S Enterococcus species composition (within-genus RA; Top", TOP_N, " + Other)"),
  ylab  = "Relative abundance within Enterococcus genus",
  stub  = paste0("stackedbar_16S_Enterococcus_species_top", TOP_N, "_other")
)

plot_stacked_bundle(
  ra_v_top,
  sample_order, sample_to_group,
  title = paste0("Enterococcus-host predicted phage composition (Host_species level; Top", TOP_N, " + Other)"),
  ylab  = "Relative abundance within Enterococcus-host phage subset",
  stub  = paste0("stackedbar_Virome_EnterococcusHostSpeciesPhage_top", TOP_N, "_other")
)

# -------------------------
# 6) Save params
# -------------------------
params <- list(
  inputs = list(
    metadata = META_FILE,
    table_16S = TABLE_16S,
    tpm = TPM_FILE,
    cherry = CHERRY_FILE
  ),
  group_order = GROUP_ORDER,
  group_colors = GROUP_COLORS,
  top_n = TOP_N,
  strip_height_ratio = STRIP_HEIGHT_RATIO,
  legend_ncol = LEGEND_NCOL
)
write_json(params, file.path(OUTDIR, "00_manifest", "params.json"), pretty = TRUE, auto_unbox = TRUE)

cat("Done. Output folder:", OUTDIR, "\n")