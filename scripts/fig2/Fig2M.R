#!/usr/bin/env Rscript

# ============================================================
# 3-group multi-omics targeted Sankey workflow from raw data
# - 16S raw table: featuretable_no_singletons.txt
# - metadata:      Sample_Metadata_3groups.csv
# - metabolomics:  Q2.xlsx (sheet = "sample")
#
# Groups:
#   A = GR
#   B = MR
#   C = MR+Phage cocktail
#
# Main outputs:
#   01_input_checks/
#   02_species_differential/
#   03_metabolite_differential/
#   04_species_metabolite_correlations/
#   05_sankey/
#
# Notes:
#   1) This script starts from raw 16S + raw metabolomics values.
#   2) 16S differential analysis uses species-level aggregation + CLR + Wilcoxon.
#   3) Metabolomics differential analysis follows your platform-style
#      0 -> half-min-positive imputation + log2 + one-way ANOVA + pooled-LSD.
#   4) Sankey is built from selected differential species/metabolites and their
#      pairwise Spearman correlations across the 15 shared samples.
#   5) Current plotting rules:
#      - metabolite side = intersection of differential metabolites from
#        GR vs MR and MR+Phage cocktail vs MR
#      - species side    = intersection of differential species from
#        GR vs MR and MR+Phage cocktail vs MR
#      - all plotted species = all intersection species with >=1 retained edge
#      - flow hue        = positive / negative correlation
#      - flow depth      = |rho|
#      - flow width      = -log10(FDR), rescaled for plotting
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggalluvial)
})

# -----------------------------
# user parameters
# -----------------------------
FEATURE_TABLE_FILE <- "featuretable_no_singletons.txt"
METADATA_FILE      <- "Sample_Metadata_3groups.csv"
METABOLOMICS_FILE  <- "Q2.xlsx"
METABOLOMICS_SHEET <- "sample"

OUTDIR <- "results_3group_species_metabolite_sankey_from_raw"

# 16S processing
SPECIES_PREVALENCE_MIN <- 0.10
SPECIES_MEAN_REL_MIN   <- 1e-04
SPECIES_CLR_PSEUDOCOUNT <- 0.5

# differential thresholds for node candidate selection
SPECIES_NODE_MODE    <- "intersect_ab_cb"   # fixed per current request
SPECIES_AB_FDR_CUT   <- 0.10
SPECIES_CB_FDR_CUT   <- 0.10

METAB_NODE_MODE      <- "intersect_ab_cb"   # fixed per current request
METAB_AB_FDR_CUT     <- 0.20
METAB_CB_FDR_CUT     <- 0.20

# all species with >=1 retained significant edge will be kept for final Sankey plot
# (multi-omics ranking is still computed, but no top-N truncation is applied)
SPECIES_TOP_N_FOR_PLOT <- NA

# correlation edge thresholds
EDGE_P_CUT           <- 0.05
EDGE_FDR_CUT         <- 0.20
EDGE_MIN_ABS_RHO     <- 0.30
EDGE_MAX_N           <- Inf   # keep all by default

# plotting
GROUP_COLORS <- c(
  "GR" = "#7CAE00",
  "MR" = "#F8766D",
  "MR+Phage cocktail" = "#00BFC4"
)

POS_LIGHT <- "#DEEBF7"
POS_DARK  <- "#4575B4"
NEG_LIGHT <- "#FEE5D9"
NEG_DARK  <- "#D73027"

# -----------------------------
# helpers
# -----------------------------
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "01_input_checks"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "02_species_differential"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "03_metabolite_differential"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "04_species_metabolite_correlations"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "05_sankey"), showWarnings = FALSE, recursive = TRUE)

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n")

stop_if_missing_file <- function(x) {
  if (!file.exists(x)) stop("File not found: ", x)
}

safe_write_csv <- function(x, file) {
  readr::write_csv(x, file, na = "")
}

clean_taxon_label <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "[^A-Za-z0-9_\\-\\+\\.]", "_")
  x <- str_replace_all(x, "_+", "_")
  x <- str_replace_all(x, "^_|_$", "")
  x
}

parse_tax_rank <- function(tax_vec, prefix) {
  out <- vapply(strsplit(as.character(tax_vec), ";\\s*"), function(parts) {
    hit <- parts[startsWith(parts, prefix)]
    if (length(hit) == 0) return(NA_character_)
    sub(paste0("^", prefix), "", hit[1])
  }, character(1))
  out[out %in% c("", NA)] <- NA_character_
  out
}

make_species_label <- function(species, genus) {
  raw <- ifelse(!is.na(species) &
                  !species %in% c(
                    "uncultured_bacterium", "unclassified_bacterium",
                    "uncultured", "unclassified", "metagenome"
                  ),
                species,
                ifelse(!is.na(genus) &
                         !genus %in% c("uncultured", "unclassified", "Ambiguous_taxa"),
                       paste0("unclassified_", genus),
                       "unclassified_species"))
  clean_taxon_label(raw)
}

clr_matrix <- function(count_mat, pseudocount = 0.5) {
  x <- as.matrix(count_mat)
  x <- x + pseudocount
  logx <- log(x)
  sweep(logx, 2, colMeans(logx), "-")
}

wilcox_p <- function(x, g, g1, g2) {
  idx <- g %in% c(g1, g2) & !is.na(x)
  if (sum(idx) < 4) return(NA_real_)
  x1 <- x[idx & g == g1]
  x2 <- x[idx & g == g2]
  if (length(x1) < 2 || length(x2) < 2) return(NA_real_)
  stats::wilcox.test(x1, x2, exact = FALSE)$p.value
}

kw_p <- function(x, g) {
  idx <- !is.na(x) & !is.na(g)
  if (sum(idx) < 6 || length(unique(g[idx])) < 2) return(NA_real_)
  stats::kruskal.test(x[idx] ~ g[idx])$p.value
}

lsd_p <- function(mean1, mean2, n1, n2, mse, df_res) {
  if (any(is.na(c(mean1, mean2, n1, n2, mse, df_res))) ||
      n1 < 2 || n2 < 2 || mse <= 0 || df_res <= 0) {
    return(NA_real_)
  }
  tval <- (mean1 - mean2) / sqrt(mse * (1 / n1 + 1 / n2))
  2 * stats::pt(abs(tval), df = df_res, lower.tail = FALSE)
}

calc_rescue_score <- function(mean_A, mean_B, mean_C) {
  denom <- abs(mean_A - mean_B)
  if (is.na(denom) || denom == 0) return(NA_real_)
  (abs(mean_A - mean_B) - abs(mean_A - mean_C)) / denom
}

choose_nodes <- function(tbl, mode = "ab_signature", ab_col, cb_col, ab_cut, cb_cut,
                         effect_ab_col = NULL, effect_cb_col = NULL, rescue_col = NULL) {
  if (mode == "ab_signature") {
    tbl %>% filter(.data[[ab_col]] < ab_cut)
  } else if (mode == "union_ab_cb") {
    tbl %>% filter((.data[[ab_col]] < ab_cut) | (.data[[cb_col]] < cb_cut))
  } else if (mode == "intersect_ab_cb") {
    tbl %>% filter((.data[[ab_col]] < ab_cut) & (.data[[cb_col]] < cb_cut))
  } else if (mode == "intersect_ab_cb_rescue") {
    if (is.null(effect_ab_col) || is.null(effect_cb_col) || is.null(rescue_col)) {
      stop("intersect_ab_cb_rescue mode requires effect_ab_col, effect_cb_col, and rescue_col")
    }
    tbl %>%
      filter(.data[[ab_col]] < ab_cut) %>%
      filter(.data[[cb_col]] < cb_cut) %>%
      filter(!is.na(.data[[effect_ab_col]]), !is.na(.data[[effect_cb_col]])) %>%
      filter(sign(.data[[effect_ab_col]]) == sign(.data[[effect_cb_col]])) %>%
      filter(!is.na(.data[[rescue_col]]), .data[[rescue_col]] > 0)
  } else {
    stop("Unsupported node mode: ", mode)
  }
}

safe_neglog10 <- function(x, floor = 1e-300) {
  -log10(pmax(x, floor, na.rm = TRUE))
}

rescale_range <- function(x, to = c(0, 1), from = NULL) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(mean(to), length(x)))
  if (is.null(from)) from <- range(x, na.rm = TRUE)
  if (!all(is.finite(from)) || diff(from) == 0) return(rep(mean(to), length(x)))
  scaled <- (x - from[1]) / diff(from)
  to[1] + scaled * diff(to)
}

interp_color <- function(col_light, col_dark, t) {
  t <- max(0, min(1, t))
  rgb1 <- grDevices::col2rgb(col_light)[, 1]
  rgb2 <- grDevices::col2rgb(col_dark)[, 1]
  rgb_new <- round(rgb1 + (rgb2 - rgb1) * t)
  grDevices::rgb(rgb_new[1], rgb_new[2], rgb_new[3], maxColorValue = 255)
}

make_manual_legend_layers <- function(y_total, neglog10_range) {
  x0 <- 2.36
  bar_w <- 0.025

  # ---------- 1) Spearman's r 色条 ----------
  bar_top <- y_total * 0.95
  bar_bot <- y_total * 0.84
  n_bar <- 80

  r_vals <- seq(-1, 1, length.out = n_bar)
  y_seq <- seq(bar_bot, bar_top, length.out = n_bar + 1)

  # 与主图保持一致：positive = blue, negative = red
  bar_cols <- vapply(r_vals, function(v) {
    if (v >= 0) {
      interp_color(POS_LIGHT, POS_DARK, abs(v))
    } else {
      interp_color(NEG_LIGHT, NEG_DARK, abs(v))
    }
  }, character(1))

  layers <- list(
    annotate("text", x = x0, y = bar_top + y_total * 0.015,
             label = "Spearman's r", hjust = 0, vjust = 0, size = 3.6)
  )

  for (i in seq_len(n_bar)) {
    layers[[length(layers) + 1]] <- annotate(
      "rect",
      xmin = x0, xmax = x0 + bar_w,
      ymin = y_seq[i], ymax = y_seq[i + 1],
      fill = bar_cols[i], color = bar_cols[i]
    )
  }

  tick_vals <- c(1.0, 0.5, 0.0, -0.5, -1.0)
  tick_pos <- bar_bot + (tick_vals - (-1)) / 2 * (bar_top - bar_bot)

  for (i in seq_along(tick_vals)) {
    layers[[length(layers) + 1]] <- annotate(
      "text",
      x = x0 + bar_w + 0.015,
      y = tick_pos[i],
      label = sprintf("%.1f", tick_vals[i]),
      hjust = 0, vjust = 0.5, size = 3.0
    )
  }

  # ---------- 2) Edge significance ----------
  sig_title_y <- bar_bot - y_total * 0.04
  line_x0 <- x0
  line_x1 <- x0 + 0.05

  q_demo <- c(0.15, 0.075, 0.025)
  neglog10_demo <- -log10(q_demo)

  width_demo <- rescale_range(
    neglog10_demo,
    to = c(1.2, 4.5),
    from = neglog10_range
  )
  width_demo <- as.numeric(width_demo)
  width_demo[!is.finite(width_demo)] <- 2.5
  width_demo <- pmax(width_demo, 0.8)

  sig_y <- c(
    sig_title_y - y_total * 0.020,
    sig_title_y - y_total * 0.040,
    sig_title_y - y_total * 0.060
  )

  layers[[length(layers) + 1]] <- annotate(
    "text",
    x = x0, y = sig_title_y,
    label = "Edge significance",
    hjust = 0, vjust = 0, size = 3.6
  )

  sig_labels <- c(
    "p < 0.05 & q >= 0.1",
    "p < 0.05 & 0.05 <= q < 0.1",
    "p < 0.05 & q < 0.05"
  )

  for (i in 1:3) {
    layers[[length(layers) + 1]] <- annotate(
      "segment",
      x = line_x0, xend = line_x1,
      y = sig_y[i], yend = sig_y[i],
      colour = "grey35",
      linewidth = width_demo[i]
    )
    layers[[length(layers) + 1]] <- annotate(
      "text",
      x = line_x1 + 0.015,
      y = sig_y[i],
      label = sig_labels[i],
      hjust = 0, vjust = 0.5, size = 3.0
    )
  }

  layers
}

write_selection_rules <- function(file,
                                  species_intersection_n,
                                  species_top_n,
                                  metab_intersection_n,
                                  plot_metab_n,
                                  all_sig_edge_n,
                                  plot_sig_edge_n) {
  txt <- c(
    "3-group species-metabolite Sankey feature selection rules",
    "========================================================",
    "",
    "Groups used:",
    "  - GR",
    "  - MR",
    "  - MR+Phage cocktail",
    "",
    "Metabolite node rules:",
    sprintf("  1) Start from all metabolites in %s (sheet = %s).", METABOLOMICS_FILE, METABOLOMICS_SHEET),
    "  2) Differential analysis uses 0 -> half-min-positive imputation, log2 transform, one-way ANOVA, and pooled-LSD pairwise comparisons.",
    sprintf("  3) Differential metabolite intersection requires BOTH (GR vs MR, FDR < %.3f) AND (MR+Phage cocktail vs MR, FDR < %.3f).",
            METAB_AB_FDR_CUT, METAB_CB_FDR_CUT),
    sprintf("  4) Number of metabolite intersection candidates = %d.", metab_intersection_n),
    sprintf("  5) Number of metabolites actually appearing in the final Sankey = %d (metabolites with at least one retained edge after keeping all plotted species).",
            plot_metab_n),
    "",
    "Species node rules:",
    "  1) 16S starts from the raw feature table and is aggregated at species level.",
    sprintf("  2) Species prefilter: prevalence >= %.2f and mean relative abundance >= %.6f.",
            SPECIES_PREVALENCE_MIN, SPECIES_MEAN_REL_MIN),
    "  3) Species differential analysis uses CLR-transformed species abundance and Wilcoxon rank-sum tests.",
    sprintf("  4) Differential species intersection requires BOTH (GR vs MR, FDR < %.3f) AND (MR+Phage cocktail vs MR, FDR < %.3f).",
            SPECIES_AB_FDR_CUT, SPECIES_CB_FDR_CUT),
    sprintf("  5) Number of species intersection candidates = %d.", species_intersection_n),
    "",
    "Multi-omics explanatory power ranking for species:",
    "  Ranking is computed within the differential-species intersection after correlating all intersection species against all differential-metabolite intersection candidates.",
    sprintf("  Correlation edge filter: raw P < %.3f, FDR < %.3f, and |rho| >= %.3f.",
            EDGE_P_CUT, EDGE_FDR_CUT, EDGE_MIN_ABS_RHO),
    "  Species-level summary metrics:",
    "    - n_linked_metabolites: number of retained significant metabolite links",
    "    - mean_abs_rho_sig: mean |rho| across retained links",
    "    - mean_neglog10_fdr_sig: mean -log10(FDR) across retained links",
    "    - species_diff_evidence: max[-log10(FDR_GR_vs_MR), -log10(FDR_Phage_vs_MR)]",
    "  Composite ranking score:",
    "    multiomics_explanatory_score = n_linked_metabolites * mean_abs_rho_sig * mean_neglog10_fdr_sig",
    "  Tie-breakers:",
    "    1) larger n_linked_metabolites",
    "    2) larger mean_abs_rho_sig",
    "    3) larger mean_neglog10_fdr_sig",
    "    4) larger species_diff_evidence",
    "  Final plotted species = all differential-intersection species with at least one retained significant edge.",
    "",
    "Flow rules in the final Sankey:",
    "  - Hue = correlation sign (positive / negative)",
    "  - Color depth = |rho|",
    "  - Flow width = -log10(FDR), rescaled for plotting",
    "",
    sprintf("Retained significant edges among all intersection species x intersection metabolites = %d.", all_sig_edge_n),
    sprintf("Retained significant edges actually plotted after keeping all plotted species = %d.", plot_sig_edge_n),
    "",
    "Output tables:",
    "  - All differential species table",
    "  - Differential-species intersection candidate table",
    "  - Species multi-omics ranking table",
    "  - All plotted species with >=1 retained edge",
    "  - All differential metabolite table",
    "  - Differential-metabolite intersection candidate table",
    "  - All correlation results for intersection species x intersection metabolites",
    "  - All retained significant edges before final plotting",
    "  - Final Sankey nodes and links used for plotting"
  )
  writeLines(txt, con = file)
}

# -----------------------------
# 1. read metadata
# -----------------------------
msg("Reading metadata...")
stop_if_missing_file(METADATA_FILE)
meta <- readr::read_csv(METADATA_FILE, show_col_types = FALSE)

required_meta_cols <- c("sample", "group")
if (!all(required_meta_cols %in% colnames(meta))) {
  stop("Metadata must contain columns: sample, group")
}

group_levels <- c("GR", "MR", "MR+Phage cocktail")
meta <- meta %>%
  filter(group %in% group_levels) %>%
  mutate(group = factor(group, levels = group_levels)) %>%
  arrange(group, sample)

samples_keep <- meta$sample

safe_write_csv(meta, file.path(OUTDIR, "01_input_checks", "metadata_used.csv"))
msg("Metadata samples: ", length(samples_keep))

# -----------------------------
# 2. read 16S feature table
# -----------------------------
msg("Reading 16S feature table...")
stop_if_missing_file(FEATURE_TABLE_FILE)
ft <- readr::read_tsv(FEATURE_TABLE_FILE, skip = 1, show_col_types = FALSE)

required_ft_cols <- c("#FeatureID", samples_keep, "taxonomy")
if (!all(required_ft_cols %in% colnames(ft))) {
  stop("16S feature table is missing required columns.")
}

ft_use <- ft %>%
  select(`#FeatureID`, all_of(samples_keep), taxonomy)

ft_tax <- ft_use %>%
  mutate(
    species = parse_tax_rank(taxonomy, "s__"),
    genus   = parse_tax_rank(taxonomy, "g__"),
    species_label = make_species_label(species, genus)
  )

safe_write_csv(
  ft_tax %>% select(`#FeatureID`, taxonomy, genus, species, species_label),
  file.path(OUTDIR, "01_input_checks", "feature_taxonomy_parsed.csv")
)

species_count <- ft_tax %>%
  select(species_label, all_of(samples_keep)) %>%
  group_by(species_label) %>%
  summarise(across(all_of(samples_keep), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

species_count_mat <- species_count %>%
  tibble::column_to_rownames("species_label") %>%
  as.matrix()

species_relab_mat <- sweep(species_count_mat, 2, colSums(species_count_mat), "/")
species_prev <- rowMeans(species_relab_mat > 0)
species_mean_relab <- rowMeans(species_relab_mat)

species_filter_tbl <- tibble::tibble(
  species = rownames(species_relab_mat),
  prevalence = species_prev,
  mean_relab = species_mean_relab,
  keep = prevalence >= SPECIES_PREVALENCE_MIN & mean_relab >= SPECIES_MEAN_REL_MIN
)

safe_write_csv(species_filter_tbl, file.path(OUTDIR, "02_species_differential", "species_filter_table.csv"))

species_keep <- species_filter_tbl %>% filter(keep) %>% pull(species)
species_count_filt <- species_count_mat[species_keep, , drop = FALSE]
species_relab_filt <- species_relab_mat[species_keep, , drop = FALSE]
species_clr_filt <- clr_matrix(species_count_filt, pseudocount = SPECIES_CLR_PSEUDOCOUNT)

species_relab_by_sample <- as.data.frame(species_relab_filt) %>%
  tibble::rownames_to_column("species")
species_clr_by_sample <- as.data.frame(species_clr_filt) %>%
  tibble::rownames_to_column("species")

safe_write_csv(species_relab_by_sample, file.path(OUTDIR, "02_species_differential", "species_relative_abundance_by_sample.csv"))
safe_write_csv(species_clr_by_sample, file.path(OUTDIR, "02_species_differential", "species_clr_by_sample.csv"))

msg("Species retained after filtering: ", nrow(species_relab_by_sample))

# -----------------------------
# 3. differential species (3 groups, CLR + Wilcoxon)
# -----------------------------
msg("Running differential species analysis...")

species_long <- as.data.frame(species_clr_filt) %>%
  tibble::rownames_to_column("species") %>%
  pivot_longer(cols = all_of(samples_keep), names_to = "sample", values_to = "clr") %>%
  left_join(meta, by = "sample")

species_diff <- species_long %>%
  group_by(species) %>%
  group_modify(~{
    d <- .x

    sum_tbl <- d %>%
      group_by(group) %>%
      summarise(
        n = sum(!is.na(clr)),
        mean = mean(clr, na.rm = TRUE),
        sd = sd(clr, na.rm = TRUE),
        median = median(clr, na.rm = TRUE),
        .groups = "drop"
      )

    getv <- function(g, field) {
      x <- sum_tbl %>% filter(group == g) %>% pull({{field}})
      if (length(x) == 0) NA_real_ else x[1]
    }

    mean_A <- getv("GR", mean)
    mean_B <- getv("MR", mean)
    mean_C <- getv("MR+Phage cocktail", mean)
    sd_A <- getv("GR", sd)
    sd_B <- getv("MR", sd)
    sd_C <- getv("MR+Phage cocktail", sd)
    med_A <- getv("GR", median)
    med_B <- getv("MR", median)
    med_C <- getv("MR+Phage cocktail", median)

    p_overall <- kw_p(d$clr, d$group)
    p_A_vs_B <- wilcox_p(d$clr, d$group, "GR", "MR")
    p_C_vs_B <- wilcox_p(d$clr, d$group, "MR+Phage cocktail", "MR")
    p_C_vs_A <- wilcox_p(d$clr, d$group, "MR+Phage cocktail", "GR")

    log2FC_A_vs_B <- mean_A - mean_B
    log2FC_C_vs_B <- mean_C - mean_B
    log2FC_C_vs_A <- mean_C - mean_A

    rescue_score <- calc_rescue_score(mean_A, mean_B, mean_C)
    direction_AB <- ifelse(is.na(log2FC_A_vs_B), NA_character_,
                           ifelse(log2FC_A_vs_B > 0, "Higher in GR", "Higher in MR"))

    tibble::tibble(
      p_overall = p_overall,
      mean_A = mean_A, sd_A = sd_A, median_A = med_A,
      mean_B = mean_B, sd_B = sd_B, median_B = med_B,
      mean_C = mean_C, sd_C = sd_C, median_C = med_C,
      log2FC_A_vs_B = log2FC_A_vs_B, p_A_vs_B = p_A_vs_B,
      log2FC_C_vs_B = log2FC_C_vs_B, p_C_vs_B = p_C_vs_B,
      log2FC_C_vs_A = log2FC_C_vs_A, p_C_vs_A = p_C_vs_A,
      rescue_score_C_toward_A = rescue_score,
      direction_AB = direction_AB
    )
  }) %>%
  ungroup() %>%
  mutate(
    fdr_overall = p.adjust(p_overall, method = "BH"),
    fdr_A_vs_B = p.adjust(p_A_vs_B, method = "BH"),
    fdr_C_vs_B = p.adjust(p_C_vs_B, method = "BH"),
    fdr_C_vs_A = p.adjust(p_C_vs_A, method = "BH")
  ) %>%
  arrange(fdr_A_vs_B, p_A_vs_B)

safe_write_csv(species_diff, file.path(OUTDIR, "02_species_differential", "species_differential_all.csv"))

species_nodes <- choose_nodes(
  tbl = species_diff,
  mode = SPECIES_NODE_MODE,
  ab_col = "fdr_A_vs_B",
  cb_col = "fdr_C_vs_B",
  ab_cut = SPECIES_AB_FDR_CUT,
  cb_cut = SPECIES_CB_FDR_CUT
) %>%
  arrange(fdr_A_vs_B, fdr_C_vs_B, species)

safe_write_csv(species_nodes, file.path(OUTDIR, "02_species_differential", "species_intersection_candidates.csv"))
safe_write_csv(species_nodes, file.path(OUTDIR, "02_species_differential", "species_nodes_selected.csv"))

msg("Selected species intersection candidates: ", nrow(species_nodes))

# -----------------------------
# 4. metabolomics raw values (3 groups only)
# -----------------------------
msg("Reading metabolomics raw table...")
stop_if_missing_file(METABOLOMICS_FILE)
met_raw <- readxl::read_excel(METABOLOMICS_FILE, sheet = METABOLOMICS_SHEET)
met_raw <- as.data.frame(met_raw)

sample_cols_met <- samples_keep
required_met_cols <- c("Compound Name (nmol/g)", sample_cols_met)
if (!all(required_met_cols %in% colnames(met_raw))) {
  stop("Metabolomics table is missing required columns in sheet '", METABOLOMICS_SHEET, "'.")
}

meta_cols_met <- intersect(
  c("中文名称", "CAS", "Formula", "Pubchem.ID", "HMDB", "KEGG", "Class", "分类", "Average Molecular Weight"),
  colnames(met_raw)
)

met_long <- met_raw %>%
  select(all_of(c("Compound Name (nmol/g)", meta_cols_met, sample_cols_met))) %>%
  pivot_longer(cols = all_of(sample_cols_met), names_to = "sample", values_to = "value_raw") %>%
  mutate(
    metabolite = as.character(`Compound Name (nmol/g)`),
    value_raw = suppressWarnings(as.numeric(value_raw))
  ) %>%
  left_join(meta, by = "sample") %>%
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

met_log2_by_sample <- met_long2 %>%
  select(metabolite, sample, value_log2) %>%
  pivot_wider(names_from = sample, values_from = value_log2)

safe_write_csv(met_log2_by_sample, file.path(OUTDIR, "03_metabolite_differential", "metabolite_log2_by_sample.csv"))

# differential metabolites
msg("Running differential metabolite analysis...")

met_diff_core <- met_long2 %>%
  group_by(metabolite) %>%
  group_modify(~{
    d <- .x

    sum_tbl <- d %>%
      group_by(group) %>%
      summarise(
        n = sum(!is.na(value_log2)),
        mean = mean(value_log2, na.rm = TRUE),
        sd = sd(value_log2, na.rm = TRUE),
        .groups = "drop"
      )

    getv <- function(g, field) {
      x <- sum_tbl %>% filter(group == g) %>% pull({{field}})
      if (length(x) == 0) NA_real_ else x[1]
    }

    mean_A <- getv("GR", mean); sd_A <- getv("GR", sd); nA <- getv("GR", n)
    mean_B <- getv("MR", mean); sd_B <- getv("MR", sd); nB <- getv("MR", n)
    mean_C <- getv("MR+Phage cocktail", mean); sd_C <- getv("MR+Phage cocktail", sd); nC <- getv("MR+Phage cocktail", n)

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

    rescue_score <- calc_rescue_score(mean_A, mean_B, mean_C)
    direction_AB <- ifelse(is.na(log2FC_A_vs_B), NA_character_,
                           ifelse(log2FC_A_vs_B > 0, "Higher in GR", "Higher in MR"))

    tibble::tibble(
      p_overall = p_overall,
      mean_A = mean_A, sd_A = sd_A,
      mean_B = mean_B, sd_B = sd_B,
      mean_C = mean_C, sd_C = sd_C,
      log2FC_A_vs_B = log2FC_A_vs_B, p_A_vs_B = p_A_vs_B,
      log2FC_C_vs_B = log2FC_C_vs_B, p_C_vs_B = p_C_vs_B,
      log2FC_C_vs_A = log2FC_C_vs_A, p_C_vs_A = p_C_vs_A,
      gm_A = gm_A, gm_B = gm_B, gm_C = gm_C,
      rescue_score_C_toward_A = rescue_score,
      direction_AB = direction_AB
    )
  }) %>%
  ungroup()

met_meta_one <- met_raw %>%
  transmute(
    metabolite = as.character(`Compound Name (nmol/g)`),
    中文名称 = if ("中文名称" %in% colnames(met_raw)) as.character(`中文名称`) else NA_character_,
    CAS = if ("CAS" %in% colnames(met_raw)) as.character(CAS) else NA_character_,
    Formula = if ("Formula" %in% colnames(met_raw)) as.character(Formula) else NA_character_,
    Pubchem.ID = if ("Pubchem.ID" %in% colnames(met_raw)) as.character(`Pubchem.ID`) else NA_character_,
    HMDB = if ("HMDB" %in% colnames(met_raw)) as.character(HMDB) else NA_character_,
    KEGG = if ("KEGG" %in% colnames(met_raw)) as.character(KEGG) else NA_character_,
    Class = if ("Class" %in% colnames(met_raw)) as.character(Class) else NA_character_,
    分类 = if ("分类" %in% colnames(met_raw)) as.character(`分类`) else NA_character_,
    `Average Molecular Weight` = if ("Average Molecular Weight" %in% colnames(met_raw)) suppressWarnings(as.numeric(`Average Molecular Weight`)) else NA_real_
  ) %>%
  distinct()

met_diff <- met_meta_one %>%
  left_join(met_diff_core, by = "metabolite") %>%
  mutate(
    fdr_overall = p.adjust(p_overall, method = "BH"),
    fdr_A_vs_B = p.adjust(p_A_vs_B, method = "BH"),
    fdr_C_vs_B = p.adjust(p_C_vs_B, method = "BH"),
    fdr_C_vs_A = p.adjust(p_C_vs_A, method = "BH")
  ) %>%
  arrange(fdr_A_vs_B, p_A_vs_B)

safe_write_csv(met_diff, file.path(OUTDIR, "03_metabolite_differential", "metabolite_differential_all.csv"))

met_nodes <- choose_nodes(
  tbl = met_diff,
  mode = METAB_NODE_MODE,
  ab_col = "fdr_A_vs_B",
  cb_col = "fdr_C_vs_B",
  ab_cut = METAB_AB_FDR_CUT,
  cb_cut = METAB_CB_FDR_CUT
) %>%
  arrange(fdr_A_vs_B, fdr_C_vs_B, metabolite)

safe_write_csv(met_nodes, file.path(OUTDIR, "03_metabolite_differential", "metabolite_intersection_candidates.csv"))
safe_write_csv(met_nodes, file.path(OUTDIR, "03_metabolite_differential", "metabolite_nodes_selected.csv"))
msg("Selected metabolite intersection candidates: ", nrow(met_nodes))

# -----------------------------
# 5. species-metabolite correlations for ALL intersection candidates
# -----------------------------
msg("Running species-metabolite correlations...")

cand_species <- intersect(species_nodes$species, rownames(species_clr_filt))
cand_metabs  <- intersect(met_nodes$metabolite, met_log2_by_sample$metabolite)

if (length(cand_species) == 0) stop("No species intersection candidates selected. Please relax species thresholds.")
if (length(cand_metabs) == 0) stop("No metabolite intersection candidates selected. Please relax metabolite thresholds.")

species_mat_cand <- species_clr_filt[cand_species, samples_keep, drop = FALSE]
met_mat_cand <- met_log2_by_sample %>%
  filter(metabolite %in% cand_metabs) %>%
  tibble::column_to_rownames("metabolite") %>%
  as.matrix()
met_mat_cand <- met_mat_cand[, samples_keep, drop = FALSE]

safe_write_csv(
  as.data.frame(species_mat_cand) %>% tibble::rownames_to_column("species"),
  file.path(OUTDIR, "04_species_metabolite_correlations", "intersection_candidate_species_clr_matrix.csv")
)
safe_write_csv(
  as.data.frame(met_mat_cand) %>% tibble::rownames_to_column("metabolite"),
  file.path(OUTDIR, "04_species_metabolite_correlations", "intersection_candidate_metabolite_log2_matrix.csv")
)

cor_long <- expand.grid(
  species = rownames(species_mat_cand),
  metabolite = rownames(met_mat_cand),
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  rowwise() %>%
  mutate(
    rho = {
      x <- as.numeric(species_mat_cand[species, ])
      y <- as.numeric(met_mat_cand[metabolite, ])
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) < 5) NA_real_ else suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
    },
    p = {
      x <- as.numeric(species_mat_cand[species, ])
      y <- as.numeric(met_mat_cand[metabolite, ])
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) < 5) {
        NA_real_
      } else {
        suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)$p.value)
      }
    }
  ) %>%
  ungroup() %>%
  mutate(
    fdr = p.adjust(p, method = "BH"),
    direction = ifelse(rho >= 0, "positive", "negative"),
    abs_rho = abs(rho)
  ) %>%
  left_join(
    species_nodes %>%
      transmute(
        species,
        species_direction_AB = direction_AB,
        species_rescue_score = rescue_score_C_toward_A,
        species_fdr_A_vs_B = fdr_A_vs_B,
        species_fdr_C_vs_B = fdr_C_vs_B
      ),
    by = "species"
  ) %>%
  left_join(
    met_nodes %>%
      transmute(
        metabolite,
        metabolite_direction_AB = direction_AB,
        metabolite_rescue_score = rescue_score_C_toward_A,
        metabolite_fdr_A_vs_B = fdr_A_vs_B,
        metabolite_fdr_C_vs_B = fdr_C_vs_B,
        中文名称
      ),
    by = "metabolite"
  ) %>%
  arrange(fdr, desc(abs_rho))

safe_write_csv(cor_long, file.path(OUTDIR, "04_species_metabolite_correlations", "species_metabolite_correlations_all_intersection_candidates.csv"))
safe_write_csv(cor_long, file.path(OUTDIR, "04_species_metabolite_correlations", "species_metabolite_correlations_all.csv"))

rho_mat <- cor_long %>%
  select(species, metabolite, rho) %>%
  pivot_wider(names_from = metabolite, values_from = rho)
p_mat <- cor_long %>%
  select(species, metabolite, p) %>%
  pivot_wider(names_from = metabolite, values_from = p)
fdr_mat <- cor_long %>%
  select(species, metabolite, fdr) %>%
  pivot_wider(names_from = metabolite, values_from = fdr)

safe_write_csv(rho_mat, file.path(OUTDIR, "04_species_metabolite_correlations", "rho_matrix.csv"))
safe_write_csv(p_mat, file.path(OUTDIR, "04_species_metabolite_correlations", "p_matrix.csv"))
safe_write_csv(fdr_mat, file.path(OUTDIR, "04_species_metabolite_correlations", "fdr_matrix.csv"))

edge_sig_all <- cor_long %>%
  filter(!is.na(rho), !is.na(p), !is.na(fdr)) %>%
  filter(p < EDGE_P_CUT, fdr < EDGE_FDR_CUT, abs_rho >= EDGE_MIN_ABS_RHO) %>%
  mutate(neglog10_fdr = safe_neglog10(fdr)) %>%
  arrange(fdr, desc(abs_rho))

if (is.finite(EDGE_MAX_N) && nrow(edge_sig_all) > EDGE_MAX_N) {
  edge_sig_all <- edge_sig_all %>% slice_head(n = EDGE_MAX_N)
}

safe_write_csv(edge_sig_all, file.path(OUTDIR, "04_species_metabolite_correlations", "significant_edges_all_intersection_candidates.csv"))
safe_write_csv(edge_sig_all, file.path(OUTDIR, "04_species_metabolite_correlations", "significant_edges_for_sankey.csv"))

msg("Significant edges retained among all intersection candidates: ", nrow(edge_sig_all))

# -----------------------------
# 6. rank species by multi-omics explanatory power within intersection candidates; keep all species with >=1 retained significant edge
# -----------------------------
msg("Ranking species by multi-omics explanatory power...")

species_rank_tbl <- species_nodes %>%
  transmute(
    species,
    direction_AB,
    rescue_score_C_toward_A,
    species_diff_evidence = pmax(safe_neglog10(fdr_A_vs_B), safe_neglog10(fdr_C_vs_B), na.rm = TRUE),
    species_fdr_A_vs_B = fdr_A_vs_B,
    species_fdr_C_vs_B = fdr_C_vs_B
  ) %>%
  left_join(
    edge_sig_all %>%
      group_by(species) %>%
      summarise(
        n_sig_edges = n(),
        n_linked_metabolites = n_distinct(metabolite),
        mean_abs_rho_sig = mean(abs_rho, na.rm = TRUE),
        max_abs_rho_sig = max(abs_rho, na.rm = TRUE),
        mean_neglog10_fdr_sig = mean(neglog10_fdr, na.rm = TRUE),
        max_neglog10_fdr_sig = max(neglog10_fdr, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "species"
  ) %>%
  mutate(
    n_sig_edges = dplyr::coalesce(n_sig_edges, 0L),
    n_linked_metabolites = dplyr::coalesce(n_linked_metabolites, 0L),
    mean_abs_rho_sig = dplyr::coalesce(mean_abs_rho_sig, 0),
    max_abs_rho_sig = dplyr::coalesce(max_abs_rho_sig, 0),
    mean_neglog10_fdr_sig = dplyr::coalesce(mean_neglog10_fdr_sig, 0),
    max_neglog10_fdr_sig = dplyr::coalesce(max_neglog10_fdr_sig, 0),
    multiomics_explanatory_score =
      n_linked_metabolites * mean_abs_rho_sig * mean_neglog10_fdr_sig
  ) %>%
  arrange(
    desc(multiomics_explanatory_score),
    desc(n_linked_metabolites),
    desc(mean_abs_rho_sig),
    desc(mean_neglog10_fdr_sig),
    desc(species_diff_evidence),
    species
  ) %>%
  mutate(multiomics_rank = row_number())

safe_write_csv(species_rank_tbl, file.path(OUTDIR, "05_sankey", "species_multiomics_explanatory_ranking.csv"))

species_nodes_plot <- species_rank_tbl %>%
  filter(n_linked_metabolites > 0)

safe_write_csv(species_nodes_plot, file.path(OUTDIR, "05_sankey", "species_all_for_sankey.csv"))

msg("Species retained for final Sankey plot: ", nrow(species_nodes_plot))

# -----------------------------
# 7. final Sankey nodes / links based on all retained intersection species
# -----------------------------
msg("Building Sankey nodes and links...")

edge_sig_plot <- edge_sig_all %>%
  filter(species %in% species_nodes_plot$species) %>%
  arrange(fdr, desc(abs_rho), species, metabolite)

safe_write_csv(edge_sig_plot, file.path(OUTDIR, "05_sankey", "plot_edges_after_species_all.csv"))

met_nodes_plot <- met_nodes %>%
  filter(metabolite %in% unique(edge_sig_plot$metabolite)) %>%
  arrange(fdr_A_vs_B, fdr_C_vs_B, metabolite)

safe_write_csv(met_nodes_plot, file.path(OUTDIR, "05_sankey", "metabolites_used_in_plot.csv"))

species_nodes2 <- species_nodes_plot %>%
  transmute(
    node_name = species,
    node_label = species,
    node_type = "Species",
    direction_AB,
    rescue_score_C_toward_A,
    multiomics_rank,
    multiomics_explanatory_score
  )

met_nodes2 <- met_nodes_plot %>%
  transmute(
    node_name = metabolite,
    node_label = metabolite,
    node_type = "Metabolite",
    direction_AB,
    rescue_score_C_toward_A = rescue_score_C_toward_A,
    multiomics_rank = NA_real_,
    multiomics_explanatory_score = NA_real_
  )

nodes_tbl <- bind_rows(species_nodes2, met_nodes2) %>%
  distinct() %>%
  arrange(node_type, node_label) %>%
  mutate(node_id = row_number() - 1L)

safe_write_csv(nodes_tbl, file.path(OUTDIR, "05_sankey", "sankey_nodes.csv"))

plot_abs_rho_min <- min(edge_sig_plot$abs_rho, na.rm = TRUE)
plot_abs_rho_max <- max(edge_sig_plot$abs_rho, na.rm = TRUE)

links_tbl <- edge_sig_plot %>%
  left_join(nodes_tbl %>% select(node_id, node_name), by = c("species" = "node_name")) %>%
  rename(source = node_id) %>%
  left_join(nodes_tbl %>% select(node_id, node_name), by = c("metabolite" = "node_name")) %>%
  rename(target = node_id) %>%
  mutate(
    edge_sign = direction,
    edge_fill = make_edge_fill(direction, abs_rho, plot_abs_rho_min, plot_abs_rho_max),
    plot_weight = rescale_range(neglog10_fdr, to = c(1.0, 5.0))
  ) %>%
  select(
    source, target, species, metabolite, 中文名称, rho, abs_rho, p, fdr,
    neglog10_fdr, plot_weight, edge_fill, edge_sign,
    species_direction_AB, metabolite_direction_AB,
    species_rescue_score, metabolite_rescue_score
  ) %>%
  arrange(fdr, desc(abs_rho), species, metabolite)

safe_write_csv(links_tbl, file.path(OUTDIR, "05_sankey", "sankey_links.csv"))

# rules text
write_selection_rules(
  file = file.path(OUTDIR, "05_sankey", "feature_selection_rules.txt"),
  species_intersection_n = nrow(species_nodes),
  species_top_n = nrow(species_nodes_plot),
  metab_intersection_n = nrow(met_nodes),
  plot_metab_n = nrow(met_nodes_plot),
  all_sig_edge_n = nrow(edge_sig_all),
  plot_sig_edge_n = nrow(edge_sig_plot)
)

# -----------------------------
# 8. static Sankey plot
# -----------------------------
msg("Plotting Sankey...")

if (nrow(links_tbl) == 0) {
  msg("No significant edges passed thresholds. Sankey plot skipped.")
} else {
  species_levels <- rev(species_nodes_plot$species)
  metabolite_levels <- met_nodes_plot$metabolite

  plot_df <- links_tbl %>%
    mutate(
      Species = factor(species, levels = species_levels),
      Metabolite = factor(metabolite, levels = metabolite_levels),
      alluvium_id = paste(species, metabolite, sep = "___")
    ) %>%
    select(alluvium_id, plot_weight, edge_fill, Species, Metabolite) %>%
    pivot_longer(
      cols = c("Species", "Metabolite"),
      names_to = "x",
      values_to = "stratum"
    ) %>%
    mutate(
      x = factor(x, levels = c("Species", "Metabolite"))
    )

  y_total <- sum(links_tbl$plot_weight, na.rm = TRUE)

  neglog10_range <- range(links_tbl$neglog10_fdr, na.rm = TRUE)
  if (!all(is.finite(neglog10_range)) || diff(neglog10_range) == 0) {
    neglog10_range <- c(1, 3)
  }

  legend_layers <- make_manual_legend_layers(
    y_total = y_total,
    neglog10_range = neglog10_range
  )
  p_sankey <- ggplot(
    plot_df,
    aes(x = x, stratum = stratum, alluvium = alluvium_id, y = plot_weight)
  ) +
    ggalluvial::geom_flow(
      aes(fill = edge_fill),
      stat = "alluvium",
      width = 0.48,
      alpha = 0.95,
      knot.pos = 0.35,
      decreasing = FALSE
    ) +
    ggalluvial::geom_stratum(
      width = 0.48,
      fill = "grey96",
      color = "grey72",
      linewidth = 0.28
    ) +
    ggplot2::geom_text(
      stat = "stratum",
      aes(label = after_stat(stratum)),
      size = 3.4,
      lineheight = 0.92
    ) +
    scale_fill_identity() +
    scale_x_discrete(limits = c("Species", "Metabolite"), expand = c(0.16, 0.58)) +
    labs(
      title = "Species-metabolite Sankey (3 groups)",
      subtitle = paste0(
         "Metabolites = differential intersection; ",
         "Species = all differential-intersection species with >=1 retained significant edge ",
         "(multi-omics explanatory ranking computed but not truncated)"
       ),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "none",
      plot.margin = margin(12, 150, 12, 28)
    ) +
    coord_cartesian(clip = "off")

  for (ly in legend_layers) {
    p_sankey <- p_sankey + ly
  }

  ggsave(
    filename = file.path(OUTDIR, "05_sankey", "species_metabolite_sankey.pdf"),
    plot = p_sankey, width = 20, height = 16
  )
  ggsave(
    filename = file.path(OUTDIR, "05_sankey", "species_metabolite_sankey.png"),
    plot = p_sankey, width = 20, height = 16, dpi = 300
  )
}

# -----------------------------
# 9. optional xlsx outputs
# -----------------------------
if (requireNamespace("openxlsx", quietly = TRUE)) {
  msg("Writing optional xlsx summary workbooks...")

  wb1 <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb1, "species_diff_all")
  openxlsx::writeData(wb1, "species_diff_all", species_diff)
  openxlsx::addWorksheet(wb1, "species_intersection_cand")
  openxlsx::writeData(wb1, "species_intersection_cand", species_nodes)
  openxlsx::saveWorkbook(
    wb1,
    file.path(OUTDIR, "02_species_differential", "species_differential_summary.xlsx"),
    overwrite = TRUE
  )

  wb2 <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb2, "metabolite_diff_all")
  openxlsx::writeData(wb2, "metabolite_diff_all", met_diff)
  openxlsx::addWorksheet(wb2, "metab_intersection_cand")
  openxlsx::writeData(wb2, "metab_intersection_cand", met_nodes)
  openxlsx::saveWorkbook(
    wb2,
    file.path(OUTDIR, "03_metabolite_differential", "metabolite_differential_summary.xlsx"),
    overwrite = TRUE
  )

  wb3 <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb3, "all_corr_intersection")
  openxlsx::writeData(wb3, "all_corr_intersection", cor_long)
  openxlsx::addWorksheet(wb3, "sig_edges_intersection")
  openxlsx::writeData(wb3, "sig_edges_intersection", edge_sig_all)
  openxlsx::addWorksheet(wb3, "species_ranking")
  openxlsx::writeData(wb3, "species_ranking", species_rank_tbl)
  openxlsx::addWorksheet(wb3, "species_all_plot")
  openxlsx::writeData(wb3, "species_all_plot", species_nodes_plot)
  openxlsx::addWorksheet(wb3, "metab_used_plot")
  openxlsx::writeData(wb3, "metab_used_plot", met_nodes_plot)
  openxlsx::addWorksheet(wb3, "plot_edges")
  openxlsx::writeData(wb3, "plot_edges", edge_sig_plot)
  openxlsx::addWorksheet(wb3, "nodes")
  openxlsx::writeData(wb3, "nodes", nodes_tbl)
  openxlsx::addWorksheet(wb3, "links")
  openxlsx::writeData(wb3, "links", links_tbl)
  openxlsx::saveWorkbook(
    wb3,
    file.path(OUTDIR, "05_sankey", "sankey_inputs_and_edges.xlsx"),
    overwrite = TRUE
  )
}

msg("Done.")
msg("Output directory: ", normalizePath(OUTDIR, winslash = "/"))
