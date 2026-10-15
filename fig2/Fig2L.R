#!/usr/bin/env Rscript

# ============================================================
# 3-group 16S species Venn + outer ring bars from raw data
# - 16S raw table: featuretable_no_singletons.txt
# - metadata:      Sample_Metadata_3groups.csv
#
# Groups:
#   A = GR
#   B = MR
#   C = MR+Phage cocktail
#
# Main outputs:
#   results_3group_16S_species_venn_ringbars/
#     01_input_checks/
#     02_species_differential/
#     03_venn_sets/
#     04_plot/
#
# Statistical framework is aligned to the species-side filtering/differential
# logic used in run_3group_species_metabolite_sankey_from_raw_v11...:
#   - species-level aggregation from raw feature table
#   - prevalence >= 0.10
#   - mean relative abundance >= 1e-4
#   - CLR with pseudocount = 0.5
#   - Wilcoxon rank-sum tests
#   - BH-FDR
#   - Venn membership based on intersection logic across:
#       * GR vs MR
#       * MR+Phage cocktail vs MR
#
# Plot design:
#   - Center: 2-set Venn (GR vs MR ; MR+Phage cocktail vs MR)
#   - Outer ring bars:
#       outer track = GR vs MR
#       inner track = MR+Phage cocktail vs MR
#   - only regions: use venn-set colors
#   - shared region: each shared species gets its own distinct color
#   - Outward  = higher in GR / MR+Phage cocktail
#   - Inward   = higher in MR
#   - Bar length = |CLR mean difference|
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
})

# -----------------------------
# user parameters
# -----------------------------
FEATURE_TABLE_FILE <- "featuretable_no_singletons.txt"
METADATA_FILE      <- "Sample_Metadata_3groups.csv"

OUTDIR <- "results_3group_16S_species_venn_ringbars"

GROUP_LEVELS <- c("GR", "MR", "MR+Phage cocktail")

# species prefilter (same as Sankey script)
SPECIES_PREVALENCE_MIN <- 0.10
SPECIES_MEAN_REL_MIN   <- 1e-04
SPECIES_CLR_PSEUDOCOUNT <- 0.5

# venn membership thresholds (same as Sankey species side)
SPECIES_AB_FDR_CUT <- 0.10   # GR vs MR
SPECIES_CB_FDR_CUT <- 0.10   # MR+Phage cocktail vs MR

# only-region colors = venn circle colors
AONLY_COLOR <- "#7CAE00"   # GR vs MR only
BONLY_COLOR <- "#00BFC4"   # MR+Phage cocktail vs MR only

# bar length scaling
BAR_MAX_LEN_OUTER <- 0.32
BAR_MAX_LEN_INNER <- 0.24

# angular ranges (degrees) for the 3 Venn regions around the outside
# left = A only ; top = shared ; right = B only
ANGLE_RANGE_AONLY  <- c(140, 240)
ANGLE_RANGE_SHARED <- c(55, 125)
ANGLE_RANGE_BONLY  <- c(-40, 40)

# plot geometry
VENN_LEFT_X   <- -0.48
VENN_RIGHT_X  <-  0.48
VENN_CENTER_Y <-  0.00
VENN_RADIUS   <-  0.78

OUTER_TRACK_BASE <- 1.72
INNER_TRACK_BASE <- 1.48

POLY_N <- 45

# -----------------------------
# helpers
# -----------------------------
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "01_input_checks"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "02_species_differential"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "03_venn_sets"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "04_plot"), showWarnings = FALSE, recursive = TRUE)

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n")
}

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

safe_neglog10 <- function(x, floor = 1e-300) {
  -log10(pmax(x, floor, na.rm = TRUE))
}

rescale_range <- function(x, to = c(0, 1), from = NULL) {
  x <- as.numeric(x)
  if (length(x) == 0) return(numeric(0))
  if (all(is.na(x))) return(rep(mean(to), length(x)))
  if (is.null(from)) from <- range(x, na.rm = TRUE)
  if (!all(is.finite(from)) || diff(from) == 0) return(rep(mean(to), length(x)))
  scaled <- (x - from[1]) / diff(from)
  to[1] + scaled * diff(to)
}

make_distinct_colors <- function(n) {
  if (n <= 0) return(character(0))

  # Wider visual separation than a simple single-pass hue wheel:
  # 1) spread hues around the full circle
  # 2) interleave odd/even positions so adjacent legend items are less similar
  # 3) alternate luminance/chroma to further increase distinguishability
  hues <- seq(15, 375, length.out = n + 1)[1:n]
  idx  <- c(seq(1, n, by = 2), seq(2, n, by = 2))
  hues <- hues[idx]

  l_vals <- rep(c(56, 72, 63), length.out = n)
  c_vals <- rep(c(110, 90, 100), length.out = n)

  grDevices::hcl(h = hues, c = c_vals, l = l_vals, fixup = TRUE)
}

make_circle_df <- function(cx, cy, r, n = 360) {
  theta <- seq(0, 2 * pi, length.out = n)
  tibble(
    x = cx + r * cos(theta),
    y = cy + r * sin(theta)
  )
}

sector_polygon <- function(theta1_deg, theta2_deg, r_inner, r_outer, n = 40) {
  th1 <- theta1_deg * pi / 180
  th2 <- theta2_deg * pi / 180
  th_outer <- seq(th1, th2, length.out = n)
  th_inner <- seq(th2, th1, length.out = n)

  tibble(
    x = c(r_outer * cos(th_outer), r_inner * cos(th_inner)),
    y = c(r_outer * sin(th_outer), r_inner * sin(th_inner))
  )
}

assign_angles <- function(df_region, theta_start_deg, theta_end_deg) {
  n <- nrow(df_region)
  if (n == 0) return(df_region)

  if (n == 1) {
    center_deg <- mean(c(theta_start_deg, theta_end_deg))
    width_deg  <- min(abs(theta_end_deg - theta_start_deg) * 0.55, 10)
    df_region$theta_center_deg <- center_deg
    df_region$theta_min_deg <- center_deg - width_deg / 2
    df_region$theta_max_deg <- center_deg + width_deg / 2
    return(df_region)
  }

  centers <- seq(theta_start_deg, theta_end_deg, length.out = n + 2)[2:(n + 1)]
  gap <- abs(diff(range(theta_start_deg, theta_end_deg))) / (n + 1)
  width_deg <- min(gap * 0.65, 9)

  df_region$theta_center_deg <- centers
  df_region$theta_min_deg <- centers - width_deg / 2
  df_region$theta_max_deg <- centers + width_deg / 2
  df_region
}

build_ring_poly <- function(plot_tbl, track = c("outer", "inner"), max_abs_fc) {
  track <- match.arg(track)

  if (nrow(plot_tbl) == 0) return(tibble())

  base_radius <- if (track == "outer") OUTER_TRACK_BASE else INNER_TRACK_BASE
  max_len     <- if (track == "outer") BAR_MAX_LEN_OUTER else BAR_MAX_LEN_INNER
  comp_col    <- if (track == "outer") "clr_diff_A_vs_B" else "clr_diff_C_vs_B"

  df_use <- plot_tbl %>%
    filter(
      (track == "outer" & venn_region %in% c("GR_vs_MR_only", "Shared")) |
        (track == "inner" & venn_region %in% c("Phage_vs_MR_only", "Shared"))
    ) %>%
    mutate(
      comp_value = .data[[comp_col]],
      comp_abs   = abs(comp_value),
      bar_len    = ifelse(
        is.na(comp_abs) | max_abs_fc <= 0,
        0,
        comp_abs / max_abs_fc * max_len
      ),
      direction  = case_when(
        is.na(comp_value) ~ "missing",
        comp_value > 0 ~ "outward",
        comp_value < 0 ~ "inward",
        TRUE ~ "flat"
      )
    )

  poly_list <- vector("list", nrow(df_use))
  idx <- 1L

  for (i in seq_len(nrow(df_use))) {
    rr <- df_use[i, ]

    if (is.na(rr$bar_len) || rr$bar_len <= 0) next

    if (rr$direction == "outward") {
      r_in  <- base_radius
      r_out <- base_radius + rr$bar_len
    } else if (rr$direction == "inward") {
      r_in  <- base_radius - rr$bar_len
      r_out <- base_radius
    } else {
      next
    }

    poly_df <- sector_polygon(
      theta1_deg = rr$theta_min_deg,
      theta2_deg = rr$theta_max_deg,
      r_inner = r_in,
      r_outer = r_out,
      n = POLY_N
    ) %>%
      mutate(
        species = rr$species,
        region = rr$venn_region,
        comparison = ifelse(track == "outer", "GR vs MR", "MR+Phage cocktail vs MR"),
        direction = rr$direction,
        value = rr$comp_value,
        abs_value = rr$comp_abs,
        color = rr$color,
        poly_id = paste(track, rr$species, sep = "__")
      )

    poly_list[[idx]] <- poly_df
    idx <- idx + 1L
  }

  bind_rows(poly_list)
}

wrap_name <- function(x, width = 28) {
  stringr::str_wrap(x, width = width)
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

meta <- meta %>%
  filter(group %in% GROUP_LEVELS) %>%
  mutate(group = factor(group, levels = GROUP_LEVELS)) %>%
  arrange(group, sample)

samples_keep <- meta$sample
safe_write_csv(meta, file.path(OUTDIR, "01_input_checks", "metadata_used.csv"))
msg("Metadata samples: ", length(samples_keep))

# -----------------------------
# 2. read 16S feature table and build species matrix
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
  column_to_rownames("species_label") %>%
  as.matrix()

species_relab_mat <- sweep(species_count_mat, 2, colSums(species_count_mat), "/")
species_prev <- rowMeans(species_relab_mat > 0)
species_mean_relab <- rowMeans(species_relab_mat)

species_filter_tbl <- tibble(
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

safe_write_csv(
  as.data.frame(species_relab_filt) %>% rownames_to_column("species"),
  file.path(OUTDIR, "02_species_differential", "species_relative_abundance_by_sample.csv")
)
safe_write_csv(
  as.data.frame(species_clr_filt) %>% rownames_to_column("species"),
  file.path(OUTDIR, "02_species_differential", "species_clr_by_sample.csv")
)

msg("Species retained after filtering: ", nrow(species_clr_filt))

# -----------------------------
# 3. differential species (same framework as Sankey species side)
# -----------------------------
msg("Running differential species analysis...")

species_long <- as.data.frame(species_clr_filt) %>%
  rownames_to_column("species") %>%
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
      x <- sum_tbl %>% filter(group == g) %>% pull({{ field }})
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

    p_A_vs_B <- wilcox_p(d$clr, d$group, "GR", "MR")
    p_C_vs_B <- wilcox_p(d$clr, d$group, "MR+Phage cocktail", "MR")
    p_C_vs_A <- wilcox_p(d$clr, d$group, "MR+Phage cocktail", "GR")

    clr_diff_A_vs_B <- mean_A - mean_B
    clr_diff_C_vs_B <- mean_C - mean_B
    clr_diff_C_vs_A <- mean_C - mean_A

    tibble(
      mean_A = mean_A, sd_A = sd_A, median_A = med_A,
      mean_B = mean_B, sd_B = sd_B, median_B = med_B,
      mean_C = mean_C, sd_C = sd_C, median_C = med_C,
      clr_diff_A_vs_B = clr_diff_A_vs_B, p_A_vs_B = p_A_vs_B,
      clr_diff_C_vs_B = clr_diff_C_vs_B, p_C_vs_B = p_C_vs_B,
      clr_diff_C_vs_A = clr_diff_C_vs_A, p_C_vs_A = p_C_vs_A,
      direction_AB = case_when(
        is.na(clr_diff_A_vs_B) ~ NA_character_,
        clr_diff_A_vs_B > 0 ~ "Higher_in_GR",
        clr_diff_A_vs_B < 0 ~ "Higher_in_MR",
        TRUE ~ "Zero"
      ),
      direction_CB = case_when(
        is.na(clr_diff_C_vs_B) ~ NA_character_,
        clr_diff_C_vs_B > 0 ~ "Higher_in_MR+Phage cocktail",
        clr_diff_C_vs_B < 0 ~ "Higher_in_MR",
        TRUE ~ "Zero"
      )
    )
  }) %>%
  ungroup() %>%
  mutate(
    fdr_A_vs_B = p.adjust(p_A_vs_B, method = "BH"),
    fdr_C_vs_B = p.adjust(p_C_vs_B, method = "BH"),
    fdr_C_vs_A = p.adjust(p_C_vs_A, method = "BH")
  ) %>%
  mutate(
    AB_sig = !is.na(fdr_A_vs_B) & fdr_A_vs_B < SPECIES_AB_FDR_CUT,
    CB_sig = !is.na(fdr_C_vs_B) & fdr_C_vs_B < SPECIES_CB_FDR_CUT,
    venn_region = case_when(
      AB_sig & CB_sig ~ "Shared",
      AB_sig & !CB_sig ~ "GR_vs_MR_only",
      !AB_sig & CB_sig ~ "Phage_vs_MR_only",
      TRUE ~ "Not_selected"
    ),
    order_key = pmax(abs(clr_diff_A_vs_B), abs(clr_diff_C_vs_B), na.rm = TRUE)
  ) %>%
  arrange(venn_region, desc(order_key), species)

safe_write_csv(species_diff, file.path(OUTDIR, "02_species_differential", "species_differential_all.csv"))

# -----------------------------
# 4. Venn sets and color map
# -----------------------------
msg("Building Venn sets...")

venn_tbl <- species_diff %>%
  filter(venn_region != "Not_selected") %>%
  mutate(region_order = case_when(
    venn_region == "GR_vs_MR_only" ~ 1L,
    venn_region == "Shared" ~ 2L,
    venn_region == "Phage_vs_MR_only" ~ 3L,
    TRUE ~ 99L
  )) %>%
  arrange(region_order, desc(order_key), species)

# -----------------------------
# 5. assign angular layout
# -----------------------------
msg("Assigning ring layout...")

venn_A <- venn_tbl %>%
  filter(venn_region == "GR_vs_MR_only") %>%
  arrange(desc(order_key), species)

venn_S <- venn_tbl %>%
  filter(venn_region == "Shared") %>%
  arrange(desc(order_key), species)

venn_B <- venn_tbl %>%
  filter(venn_region == "Phage_vs_MR_only") %>%
  arrange(desc(order_key), species)

venn_A <- assign_angles(venn_A, ANGLE_RANGE_AONLY[1], ANGLE_RANGE_AONLY[2])
venn_S <- assign_angles(venn_S, ANGLE_RANGE_SHARED[1], ANGLE_RANGE_SHARED[2])
venn_B <- assign_angles(venn_B, ANGLE_RANGE_BONLY[1], ANGLE_RANGE_BONLY[2])

# Shared colors are assigned AFTER angular layout so that:
#   ring order (left -> right across the shared arc),
#   legend order (top -> bottom),
#   and color assignment all correspond one-to-one.
shared_species_order <- venn_S %>%
  arrange(desc(theta_center_deg)) %>%
  pull(species)

shared_color_vec <- make_distinct_colors(length(shared_species_order))
names(shared_color_vec) <- shared_species_order

color_map_tbl <- bind_rows(
  venn_A %>%
    transmute(species, venn_region, theta_center_deg, color = AONLY_COLOR),
  venn_S %>%
    transmute(species, venn_region, theta_center_deg, color = unname(shared_color_vec[species])),
  venn_B %>%
    transmute(species, venn_region, theta_center_deg, color = BONLY_COLOR)
)

venn_tbl <- venn_tbl %>%
  left_join(
    color_map_tbl %>% select(species, venn_region, theta_center_deg, color),
    by = c("species", "venn_region")
  )

venn_A <- venn_A %>%
  left_join(color_map_tbl %>% select(species, venn_region, color),
            by = c("species", "venn_region"))

venn_S <- venn_S %>%
  left_join(color_map_tbl %>% select(species, venn_region, color),
            by = c("species", "venn_region"))

venn_B <- venn_B %>%
  left_join(color_map_tbl %>% select(species, venn_region, color),
            by = c("species", "venn_region"))

safe_write_csv(venn_tbl, file.path(OUTDIR, "03_venn_sets", "species_venn_assignment_all.csv"))
safe_write_csv(color_map_tbl, file.path(OUTDIR, "03_venn_sets", "species_color_map.csv"))
safe_write_csv(
  venn_tbl %>% filter(venn_region == "GR_vs_MR_only"),
  file.path(OUTDIR, "03_venn_sets", "GR_vs_MR_only.csv")
)
safe_write_csv(
  venn_tbl %>% filter(venn_region == "Shared") %>% arrange(desc(theta_center_deg)),
  file.path(OUTDIR, "03_venn_sets", "shared.csv")
)
safe_write_csv(
  venn_tbl %>% filter(venn_region == "Phage_vs_MR_only"),
  file.path(OUTDIR, "03_venn_sets", "Phage_vs_MR_only.csv")
)

plot_layout_tbl <- bind_rows(venn_A, venn_S, venn_B) %>%
  mutate(
    label = wrap_name(species, width = 28),
    label_wrapped = label
  )

safe_write_csv(plot_layout_tbl, file.path(OUTDIR, "04_plot", "ring_bar_layout.csv"))

# -----------------------------
# 6. build ring polygons
# -----------------------------
msg("Building ring bar polygons...")

max_abs_fc <- max(
  abs(c(plot_layout_tbl$clr_diff_A_vs_B, plot_layout_tbl$clr_diff_C_vs_B)),
  na.rm = TRUE
)
if (!is.finite(max_abs_fc) || max_abs_fc <= 0) max_abs_fc <- 1

outer_poly <- build_ring_poly(plot_layout_tbl, track = "outer", max_abs_fc = max_abs_fc)
inner_poly <- build_ring_poly(plot_layout_tbl, track = "inner", max_abs_fc = max_abs_fc)
ring_poly_tbl <- bind_rows(outer_poly, inner_poly)

safe_write_csv(ring_poly_tbl, file.path(OUTDIR, "04_plot", "ring_bar_segments.csv"))

ring_base_outer <- make_circle_df(0, 0, OUTER_TRACK_BASE, n = 400)
ring_base_inner <- make_circle_df(0, 0, INNER_TRACK_BASE, n = 400)
venn_circle_left  <- make_circle_df(VENN_LEFT_X, VENN_CENTER_Y, VENN_RADIUS, n = 400)
venn_circle_right <- make_circle_df(VENN_RIGHT_X, VENN_CENTER_Y, VENN_RADIUS, n = 400)

# -----------------------------
# 7. plot
# -----------------------------
msg("Plotting Venn + outer ring bars...")

count_Aonly  <- sum(plot_layout_tbl$venn_region == "GR_vs_MR_only")
count_Shared <- sum(plot_layout_tbl$venn_region == "Shared")
count_Bonly  <- sum(plot_layout_tbl$venn_region == "Phage_vs_MR_only")

fill_values <- setNames(plot_layout_tbl$color, plot_layout_tbl$species)
legend_df <- plot_layout_tbl %>%
  filter(venn_region == "Shared") %>%
  distinct(species, label, color, theta_center_deg) %>%
  arrange(desc(theta_center_deg))
legend_breaks <- legend_df$species
legend_labels <- legend_df$label

subtitle_text <- paste0(
  "Species-level 16S; prevalence >= ", SPECIES_PREVALENCE_MIN,
  ", mean rel. abundance >= ", SPECIES_MEAN_REL_MIN,
  "; CLR pseudocount = ", SPECIES_CLR_PSEUDOCOUNT,
  "; Wilcoxon + BH-FDR. Outer track = GR vs MR; inner track = MR+Phage cocktail vs MR. ",
  "Outward = higher in the first group of each comparison; inward = higher in MR. ",
  "Bar length = |CLR mean difference|."
)

p <- ggplot() +
  geom_path(data = venn_circle_left, aes(x, y), color = "#7CAE00", linewidth = 0.8) +
  geom_path(data = venn_circle_right, aes(x, y), color = "#00BFC4", linewidth = 0.8) +
  annotate(
    "path",
    x = venn_circle_left$x, y = venn_circle_left$y,
    colour = grDevices::adjustcolor("#7CAE00", alpha.f = 0.85),
    linewidth = 0.8
  ) +
  annotate(
    "path",
    x = venn_circle_right$x, y = venn_circle_right$y,
    colour = grDevices::adjustcolor("#00BFC4", alpha.f = 0.85),
    linewidth = 0.8
  ) +
  geom_path(data = ring_base_outer, aes(x, y), color = "grey85", linewidth = 0.4) +
  geom_path(data = ring_base_inner, aes(x, y), color = "grey88", linewidth = 0.4) +
  geom_polygon(
    data = ring_poly_tbl,
    aes(x = x, y = y, group = poly_id, fill = species),
    color = NA,
    alpha = 0.95
  ) +
  annotate("text", x = VENN_LEFT_X, y = VENN_RADIUS + 0.22,
           label = "GR vs MR", color = "#4D7E00", fontface = "bold", size = 5) +
  annotate("text", x = VENN_RIGHT_X, y = VENN_RADIUS + 0.22,
           label = "MR+Phage cocktail vs MR", color = "#008CA0", fontface = "bold", size = 5) +
  annotate("text", x = -0.88, y = 0.02, label = as.character(count_Aonly), size = 7, fontface = "bold") +
  annotate("text", x =  0.00, y = 0.02, label = as.character(count_Shared), size = 7, fontface = "bold") +
  annotate("text", x =  0.88, y = 0.02, label = as.character(count_Bonly), size = 7, fontface = "bold") +
  annotate("text", x = -0.88, y = -0.22, label = "only", size = 3.8, color = "grey30") +
  annotate("text", x =  0.00, y = -0.22, label = "shared", size = 3.8, color = "grey30") +
  annotate("text", x =  0.88, y = -0.22, label = "only", size = 3.8, color = "grey30") +
  scale_fill_manual(
    values = fill_values,
    breaks = legend_breaks,
    labels = legend_labels,
    guide = guide_legend(
      ncol = 1,
      byrow = TRUE,
      keyheight = grid::unit(0.52, "cm"),
      keywidth = grid::unit(0.52, "cm")
    )
  ) +
  coord_equal(xlim = c(-2.35, 2.95), ylim = c(-2.20, 2.20), clip = "off") +
  labs(
    title = "Differential 16S species Venn with outer ring bars",
    subtitle = subtitle_text,
    fill = "Shared species"
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 8.5),
    plot.margin = margin(12, 16, 12, 16)
  )

ggsave(
  filename = file.path(OUTDIR, "04_plot", "species_venn_outer_ringbars.pdf"),
  plot = p, width = 14, height = 10
)
ggsave(
  filename = file.path(OUTDIR, "04_plot", "species_venn_outer_ringbars.png"),
  plot = p, width = 14, height = 10, dpi = 300
)

msg("Done.")
msg("Output directory: ", normalizePath(OUTDIR, winslash = "/"))