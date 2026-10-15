suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tibble)
  library(scales)
})

# ============================================================
# plot_16S_group_stacked_bars_global_and_Enterococcus_lineage.R
# Purpose:
#   1) Plot group-level stacked barplots for the whole 16S community
#      at 6 taxonomic ranks:
#         Phylum / Class / Order / Family / Genus / Species
#   2) Plot Enterococcus-lineage zoom-in stacked barplots (scheme 1):
#         - Within Enterococcus-containing Phylum -> Class composition
#         - Within Enterococcus-containing Class  -> Order composition
#         - Within Enterococcus-containing Order  -> Family composition
#         - Within Enterococcus-containing Family -> Genus composition
#         - Within Enterococcus genus            -> Species composition
# Notes:
#   - No per-sample barplots. Only group-level summaries are plotted.
#   - Relative abundance is first computed at sample level, then averaged
#     within each group.
#   - Zoom-in lineage plots are renormalized WITHIN the retained parent clade,
#     so each bar sums to 100% inside that parent taxonomic scope.
#   - Low-abundance taxa are collapsed into "Others" for cleaner figures.
# Inputs (working directory):
#   - phyloseq_16S.rds                  [required]
#   - Sample_Metadata_3groups.csv       [optional fallback]
# Outputs:
#   results_16S_group_stacked_bars_scheme1/
#     01_tables/
#     02_plots_pdf/
#     03_plots_png/
# ============================================================

# -----------------------------
# User parameters
# -----------------------------
PHYLOSEQ_RDS <- "phyloseq_16S.rds"
META_CSV     <- "Sample_Metadata_3groups.csv"
OUTDIR       <- "results_16S_group_stacked_bars_scheme1"

GROUP_LEVELS <- c("GR", "MR", "MR+Phage cocktail")
GROUP_COLORS <- c(
  "GR" = "#7CAE00",
  "MR" = "#F8766D",
  "MR+Phage cocktail" = "#00BFC4"
)

GLOBAL_RANKS <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
FOCUS_RANKS  <- c("Class", "Order", "Family", "Genus", "Species")
FOCUS_PARENTS <- c("Phylum", "Class", "Order", "Family", "Genus")

TOP_N_GLOBAL <- c(
  "Phylum"  = 10,
  "Class"   = 12,
  "Order"   = 12,
  "Family"  = 12,
  "Genus"   = 12,
  "Species" = 12
)

TOP_N_FOCUS <- c(
  "Class"   = 10,
  "Order"   = 10,
  "Family"  = 10,
  "Genus"   = 12,
  "Species" = 12
)

PDF_WIDTH  <- 7.2
PDF_HEIGHT <- 4.8
PNG_WIDTH  <- 2200
PNG_HEIGHT <- 1500
PNG_RES    <- 300

OTHERS_LABEL <- "Others"
UNCLASSIFIED_PREFIX <- "Unclassified"

# -----------------------------
# Basic checks and folders
# -----------------------------
stopifnot(file.exists(PHYLOSEQ_RDS))

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "01_tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "02_plots_pdf"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "03_plots_png"), showWarnings = FALSE, recursive = TRUE)

message("[1/7] Reading phyloseq object...")
ps <- readRDS(PHYLOSEQ_RDS)
if (!inherits(ps, "phyloseq")) {
  stop("Input file is not a phyloseq object: ", PHYLOSEQ_RDS)
}

# -----------------------------
# Helper functions
# -----------------------------
clean_tax_value <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- stringr::str_replace(x, "^[a-z]__", "")
  x[x %in% c("", "NA", "N/A", "NULL", "na", "nan")] <- NA_character_
  x[stringr::str_detect(stringr::str_to_lower(x), "unclassified|uncultured|unknown|metagenome")] <- NA_character_
  x
}

safe_filename <- function(x) {
  x %>%
    stringr::str_replace_all("[[:space:]]+", "_") %>%
    stringr::str_replace_all("[+/\\\\():;,.]", "_") %>%
    stringr::str_replace_all("_+", "_")
}

make_palette <- function(labels, others_label = OTHERS_LABEL) {
  labels <- unique(as.character(labels))
  labels_no_others <- setdiff(labels, others_label)
  n <- length(labels_no_others)

  if (n <= 0) {
    pal <- c()
  } else {
    pal <- grDevices::hcl.colors(
      n = n,
      palette = if (n <= 8) "Dark 3" else "Dynamic",
      alpha = 1
    )
    names(pal) <- labels_no_others
  }

  if (others_label %in% labels) {
    pal <- c(pal, setNames("#BDBDBD", others_label))
  }
  pal[labels]
}


phyloseq_df <- function(x) {
  if (inherits(x, "sample_data")) {
    rn <- rownames(x)
    out <- data.frame(unclass(x), stringsAsFactors = FALSE, check.names = FALSE)
    if (!is.null(rn)) rownames(out) <- rn
    return(out)
  }

  if (inherits(x, "taxonomyTable")) {
    out <- methods::as(x, "matrix")
    out <- base::as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
    return(out)
  }

  if (is.matrix(x)) {
    out <- base::as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
    return(out)
  }

  out <- base::as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  attributes(out) <- attributes(out)[intersect(names(attributes(out)), c("names", "row.names", "class"))]
  class(out) <- "data.frame"
  out
}

find_group_column <- function(meta_df) {
  nms <- colnames(meta_df)
  nms_low <- tolower(nms)
  cand <- c("group", "grp", "treatment", "condition")
  hit <- nms[match(cand, nms_low, nomatch = 0)]
  hit <- hit[hit != ""]
  if (length(hit) > 0) {
    return(hit[1])
  }
  NULL
}

standardize_group <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  out <- dplyr::case_when(
    x %in% c("GR", "Good responder", "Good responders") ~ "GR",
    x %in% c("MR", "Moderate responder", "Moderate responders") ~ "MR",
    x %in% c("MR+Phage cocktail", "MR + Phage cocktail", "MR+Phage", "Phage cocktail") ~ "MR+Phage cocktail",
    TRUE ~ x
  )
  out
}

get_sample_metadata <- function(ps, meta_csv = NULL) {
  meta_ps <- phyloseq_df(sample_data(ps))
  meta_ps$sample <- rownames(meta_ps)

  group_col <- find_group_column(meta_ps)
  if (!is.null(group_col)) {
    meta <- data.frame(
      sample = as.character(meta_ps$sample),
      group  = standardize_group(meta_ps[[group_col]]),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else if (!is.null(meta_csv) && file.exists(meta_csv)) {
    meta_csv_df <- readr::read_csv(meta_csv, show_col_types = FALSE)
    meta_csv_df <- base::as.data.frame(meta_csv_df, stringsAsFactors = FALSE, check.names = FALSE)
    nms_low <- tolower(colnames(meta_csv_df))
    sample_col <- colnames(meta_csv_df)[match("sample", nms_low, nomatch = 0)]
    group_col2 <- colnames(meta_csv_df)[match("group", nms_low, nomatch = 0)]
    if (length(sample_col) == 0 || length(group_col2) == 0) {
      stop("Metadata CSV must contain columns named sample and group.")
    }
    meta <- data.frame(
      sample = as.character(meta_csv_df[[sample_col[1]]]),
      group  = standardize_group(meta_csv_df[[group_col2[1]]]),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    stop("Cannot find group information in phyloseq sample_data, and no usable metadata CSV fallback is available.")
  }

  meta <- meta[meta$sample %in% sample_names(ps), , drop = FALSE]
  meta <- meta[meta$group %in% GROUP_LEVELS, , drop = FALSE]
  meta$group <- factor(meta$group, levels = GROUP_LEVELS)

  if (nrow(meta) == 0) {
    stop("No samples retained after matching GROUP_LEVELS. Please check group labels.")
  }

  rownames(meta) <- NULL
  meta
}

get_rank_mapping <- function(ps, wanted_ranks = GLOBAL_RANKS) {
  tt <- phyloseq_df(tax_table(ps))
  tax_cols <- colnames(tt)
  tax_cols_low <- tolower(tax_cols)

  mapped <- sapply(wanted_ranks, function(rk) {
    hit <- tax_cols[match(tolower(rk), tax_cols_low, nomatch = 0)]
    if (length(hit) == 0) NA_character_ else hit[1]
  }, USE.NAMES = TRUE)

  if (any(is.na(mapped))) {
    missing_ranks <- names(mapped)[is.na(mapped)]
    stop(
      "Missing taxonomy ranks in tax_table: ",
      paste(missing_ranks, collapse = ", "),
      "\nAvailable columns: ", paste(tax_cols, collapse = ", ")
    )
  }

  mapped
}

subset_phyloseq_samples <- function(ps, meta) {
  prune_samples(meta$sample, ps)
}

transform_to_relative <- function(ps) {
  transform_sample_counts(ps, function(x) {
    s <- sum(x, na.rm = TRUE)
    if (s <= 0) return(x)
    x / s
  })
}

ps_to_long <- function(ps, rank_col, meta, rank_label = NULL) {
  if (is.null(rank_label)) rank_label <- rank_col

  ps_rank <- tax_glom(ps, taxrank = rank_col, NArm = FALSE)

  otu <- as(otu_table(ps_rank), "matrix")
  if (!taxa_are_rows(ps_rank)) {
    otu <- t(otu)
  }

  tt <- as.data.frame(tax_table(ps_rank), stringsAsFactors = FALSE)
  tt$TaxonID <- rownames(tt)
  tt[[rank_col]] <- clean_tax_value(tt[[rank_col]])
  tt[[rank_col]][is.na(tt[[rank_col]])] <- paste0(UNCLASSIFIED_PREFIX, "_", rank_label)

  df_long <- as.data.frame(otu, stringsAsFactors = FALSE) %>%
    tibble::rownames_to_column("TaxonID") %>%
    tidyr::pivot_longer(-TaxonID, names_to = "sample", values_to = "abundance") %>%
    dplyr::left_join(tt[, c("TaxonID", rank_col), drop = FALSE], by = "TaxonID") %>%
    dplyr::left_join(meta, by = "sample")

  names(df_long)[names(df_long) == rank_col] <- "taxon"

  df_long
}

collapse_low_abundance_taxa <- function(df_group, top_n, others_label = OTHERS_LABEL) {
  overall <- df_group %>%
    dplyr::group_by(taxon) %>%
    dplyr::summarise(mean_abundance_all = mean(mean_abundance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_abundance_all), taxon)

  if (nrow(overall) > top_n) {
    keep_taxa <- overall$taxon[seq_len(top_n)]
  } else {
    keep_taxa <- overall$taxon
  }

  out <- df_group %>%
    dplyr::mutate(taxon_plot = ifelse(taxon %in% keep_taxa, taxon, others_label)) %>%
    dplyr::group_by(group, taxon_plot) %>%
    dplyr::summarise(mean_abundance = sum(mean_abundance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::rename(taxon = taxon_plot)

  ord <- out %>%
    dplyr::group_by(taxon) %>%
    dplyr::summarise(mean_abundance_all = mean(mean_abundance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_abundance_all), taxon) %>%
    dplyr::pull(taxon)

  out$taxon <- factor(out$taxon, levels = ord)
  out$group <- factor(out$group, levels = GROUP_LEVELS)
  out
}

summarize_group_mean <- function(df_long) {
  df_long %>%
    dplyr::group_by(group, taxon) %>%
    dplyr::summarise(mean_abundance = mean(abundance, na.rm = TRUE), .groups = "drop")
}

make_plot <- function(df_plot, title_text, subtitle_text = NULL, palette = NULL) {
  taxa_levels <- levels(df_plot$taxon)
  if (is.null(palette)) {
    palette <- make_palette(taxa_levels)
  }

  p <- ggplot(df_plot, aes(x = group, y = mean_abundance, fill = taxon)) +
    geom_col(width = 0.70, color = "white", linewidth = 0.35) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_fill_manual(values = palette, drop = FALSE) +
    labs(
      x = NULL,
      y = "Mean relative abundance",
      title = title_text,
      subtitle = subtitle_text,
      fill = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "#E6E6E6", linewidth = 0.35),
      axis.text.x = element_text(color = "black", size = 11),
      axis.text.y = element_text(color = "black", size = 10),
      axis.title.y = element_text(color = "black", size = 11),
      plot.title = element_text(face = "bold", size = 12, hjust = 0),
      plot.subtitle = element_text(size = 9.5, color = "#4D4D4D", hjust = 0),
      legend.position = "right",
      legend.text = element_text(size = 9),
      legend.key.height = grid::unit(0.42, "cm"),
      panel.border = element_rect(color = "black", linewidth = 0.6),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.6)
    )

  p
}

save_plot_dual <- function(plot_obj, pdf_file, png_file,
                           pdf_width = PDF_WIDTH, pdf_height = PDF_HEIGHT,
                           png_width = PNG_WIDTH, png_height = PNG_HEIGHT, png_res = PNG_RES) {
  ggsave(pdf_file, plot = plot_obj, width = pdf_width, height = pdf_height, units = "in", device = cairo_pdf)
  ggsave(png_file, plot = plot_obj, width = png_width, height = png_height, units = "px", dpi = png_res, bg = "white")
}

infer_enterococcus_lineage <- function(ps, rank_map) {
  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) {
    otu <- t(otu)
  }
  total_abund <- rowSums(otu, na.rm = TRUE)

  tt <- phyloseq_df(tax_table(ps))
  for (rk in names(rank_map)) {
    tt[[rank_map[[rk]]]] <- clean_tax_value(tt[[rank_map[[rk]]]])
  }

  idx <- which(tt[[rank_map[["Genus"]]]] == "Enterococcus")
  if (length(idx) == 0) {
    stop("Cannot find Genus == 'Enterococcus' in tax_table.")
  }

  choose_dominant <- function(values, weights) {
    keep <- !is.na(values) & values != ""
    values <- values[keep]
    weights <- weights[keep]
    if (length(values) == 0) return(NA_character_)
    agg <- tapply(weights, values, sum, na.rm = TRUE)
    names(which.max(agg))[1]
  }

  lineage <- c(
    Phylum = choose_dominant(tt[idx, rank_map[["Phylum"]]], total_abund[idx]),
    Class  = choose_dominant(tt[idx, rank_map[["Class"]]],  total_abund[idx]),
    Order  = choose_dominant(tt[idx, rank_map[["Order"]]],  total_abund[idx]),
    Family = choose_dominant(tt[idx, rank_map[["Family"]]], total_abund[idx]),
    Genus  = "Enterococcus"
  )

  as.list(lineage)
}

prune_by_parent <- function(ps, rank_col, parent_value) {
  tt <- phyloseq_df(tax_table(ps))
  vals <- clean_tax_value(tt[[rank_col]])
  keep_taxa <- rownames(tt)[which(vals == parent_value)]
  if (length(keep_taxa) == 0) {
    stop("No taxa retained for parent filter: ", rank_col, " == ", parent_value)
  }
  prune_taxa(keep_taxa, ps)
}

run_one_panel <- function(ps_rel, meta, rank_col, rank_label,
                          scope_label, top_n, out_prefix,
                          subtitle_text = NULL) {
  df_long <- ps_to_long(ps_rel, rank_col = rank_col, meta = meta, rank_label = rank_label)

  df_group <- summarize_group_mean(df_long)
  df_plot <- collapse_low_abundance_taxa(df_group, top_n = top_n, others_label = OTHERS_LABEL)

  palette <- make_palette(levels(df_plot$taxon))
  plot_title <- paste0(scope_label, " | ", rank_label)
  p <- make_plot(df_plot, title_text = plot_title, subtitle_text = subtitle_text, palette = palette)

  readr::write_csv(df_long, file.path(OUTDIR, "01_tables", paste0(out_prefix, "__per_sample_long.csv")))
  readr::write_csv(df_group, file.path(OUTDIR, "01_tables", paste0(out_prefix, "__group_mean_full.csv")))
  readr::write_csv(
    df_plot %>% dplyr::mutate(taxon = as.character(taxon)),
    file.path(OUTDIR, "01_tables", paste0(out_prefix, "__group_mean_plotted.csv"))
  )

  save_plot_dual(
    p,
    pdf_file = file.path(OUTDIR, "02_plots_pdf", paste0(out_prefix, ".pdf")),
    png_file = file.path(OUTDIR, "03_plots_png", paste0(out_prefix, ".png"))
  )

  p
}

# -----------------------------
# Prepare phyloseq + metadata
# -----------------------------
message("[2/7] Preparing metadata and checking taxonomy ranks...")
rank_map <- get_rank_mapping(ps, GLOBAL_RANKS)
meta <- get_sample_metadata(ps, meta_csv = META_CSV)
ps <- subset_phyloseq_samples(ps, meta)

sample_data(ps)$group <- factor(meta$group[match(sample_names(ps), meta$sample)], levels = GROUP_LEVELS)

meta_out <- data.frame(
  sample = sample_names(ps),
  group = as.character(sample_data(ps)$group),
  stringsAsFactors = FALSE
)
readr::write_csv(meta_out, file.path(OUTDIR, "01_tables", "samples_used_and_groups.csv"))

# -----------------------------
# Whole-community plots
# -----------------------------
message("[3/7] Building whole-community relative-abundance panels...")
ps_rel_global <- transform_to_relative(ps)

plots_global <- list()
for (rk in GLOBAL_RANKS) {
  rk_col <- rank_map[[rk]]
  out_prefix <- paste0("global__", safe_filename(rk))

  plots_global[[rk]] <- run_one_panel(
    ps_rel = ps_rel_global,
    meta = meta,
    rank_col = rk_col,
    rank_label = rk,
    scope_label = "Whole community",
    top_n = unname(TOP_N_GLOBAL[[rk]]),
    out_prefix = out_prefix,
    subtitle_text = "Sample-level relative abundance averaged within each group"
  )
}

# -----------------------------
# Infer Enterococcus lineage
# -----------------------------
message("[4/7] Inferring Enterococcus lineage...")
lineage <- infer_enterococcus_lineage(ps_rel_global, rank_map = rank_map)
lineage_df <- data.frame(
  rank = names(lineage),
  value = unlist(lineage, use.names = FALSE),
  stringsAsFactors = FALSE
)
readr::write_csv(lineage_df, file.path(OUTDIR, "01_tables", "Enterococcus_inferred_lineage.csv"))

# -----------------------------
# Focused Enterococcus-lineage plots
# -----------------------------
message("[5/7] Building Enterococcus-lineage zoom-in panels...")
plots_focus <- list()

for (i in seq_along(FOCUS_RANKS)) {
  child_rank <- FOCUS_RANKS[i]
  parent_rank <- FOCUS_PARENTS[i]

  parent_value <- unname(lineage[[parent_rank]])
  if (is.na(parent_value) || parent_value == "") {
    warning("Skipping focus plot because inferred parent is missing: ", parent_rank)
    next
  }

  parent_col <- rank_map[[parent_rank]]
  child_col  <- rank_map[[child_rank]]

  ps_parent <- prune_by_parent(ps, rank_col = parent_col, parent_value = parent_value)
  ps_parent_rel <- transform_to_relative(ps_parent)

  scope_label <- paste0("Within ", parent_rank, ": ", parent_value)
  out_prefix <- paste0(
    "focus__within_", safe_filename(parent_rank), "_", safe_filename(parent_value),
    "__", safe_filename(child_rank)
  )

  subtitle_text <- paste0(
    "Renormalized within ", parent_rank, " = ", parent_value,
    "; sample-level relative abundance averaged within each group"
  )

  plots_focus[[child_rank]] <- run_one_panel(
    ps_rel = ps_parent_rel,
    meta = meta,
    rank_col = child_col,
    rank_label = child_rank,
    scope_label = scope_label,
    top_n = unname(TOP_N_FOCUS[[child_rank]]),
    out_prefix = out_prefix,
    subtitle_text = subtitle_text
  )
}

# -----------------------------
# Multi-page summary PDFs
# -----------------------------
message("[6/7] Exporting multi-page summary PDFs...")
pdf(file.path(OUTDIR, "02_plots_pdf", "all_global_panels_multipage.pdf"), width = PDF_WIDTH, height = PDF_HEIGHT, onefile = TRUE)
for (nm in names(plots_global)) print(plots_global[[nm]])
dev.off()

pdf(file.path(OUTDIR, "02_plots_pdf", "all_Enterococcus_lineage_panels_multipage.pdf"), width = PDF_WIDTH, height = PDF_HEIGHT, onefile = TRUE)
for (nm in names(plots_focus)) print(plots_focus[[nm]])
dev.off()

# -----------------------------
# Done
# -----------------------------
message("[7/7] Done.")
message("Output directory: ", OUTDIR)
message("Tables: ", file.path(OUTDIR, "01_tables"))
message("PDF plots: ", file.path(OUTDIR, "02_plots_pdf"))
message("PNG plots: ", file.path(OUTDIR, "03_plots_png"))
message("Inferred Enterococcus lineage: ", paste(lineage_df$rank, lineage_df$value, sep = "=", collapse = " ; "))
