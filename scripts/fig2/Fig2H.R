#!/usr/bin/env Rscript

# ============================================================
# 3-group metabolomics Venn + outer ring bars from raw data
# - metadata:      Sample_Metadata_3groups.csv
# - metabolomics:  Q2.xlsx (sheet = "sample")
#
# Groups:
#   A = GR
#   B = MR
#   C = MR+Phage cocktail
#
# Main outputs:
#   results_3group_metabolite_venn_ringbars/
#     01_input_checks/
#     02_differential/
#     03_venn_sets/
#     04_plot/
#
# Plot design:
#   - Center: 2-set Venn (GR vs MR ; MR+Phage cocktail vs MR)
#   - Outer ring bars:
#       outer track = GR vs MR
#       inner track = MR+Phage cocktail vs MR
#   - only regions: use venn-set colors
#   - shared region: each shared metabolite gets its own distinct color
#   - Outward  = higher in GR / MR+Phage cocktail
#   - Inward   = lower in GR / MR+Phage cocktail
#   - Bar length = |log2FC|
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

# -----------------------------
# user parameters
# -----------------------------
METADATA_FILE      <- "Sample_Metadata_3groups.csv"
METABOLOMICS_FILE  <- "Q2.xlsx"
METABOLOMICS_SHEET <- "sample"

OUTDIR <- "results_3group_metabolite_venn_ringbars"

GROUP_LEVELS <- c("GR", "MR", "MR+Phage cocktail")

# metabolite differential cutoffs for set membership
METAB_AB_FDR_CUT <- 0.20   # GR vs MR
METAB_CB_FDR_CUT <- 0.20   # MR+Phage cocktail vs MR

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
dir.create(file.path(OUTDIR, "02_differential"), showWarnings = FALSE, recursive = TRUE)
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

lsd_p <- function(mean1, mean2, n1, n2, mse, df_res) {
  if (any(is.na(c(mean1, mean2, n1, n2, mse, df_res))) ||
      n1 < 2 || n2 < 2 || mse <= 0 || df_res <= 0) {
    return(NA_real_)
  }
  tval <- (mean1 - mean2) / sqrt(mse * (1 / n1 + 1 / n2))
  2 * stats::pt(abs(tval), df = df_res, lower.tail = FALSE)
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
  hues <- seq(15, 375, length.out = n + 1)[1:n]
  grDevices::hcl(h = hues, c = 95, l = 65)
}

make_circle_df <- function(cx, cy, r, n = 360) {
  theta <- seq(0, 2 * pi, length.out = n)
  tibble::tibble(
    x = cx + r * cos(theta),
    y = cy + r * sin(theta)
  )
}

sector_polygon <- function(theta1_deg, theta2_deg, r_inner, r_outer, n = 40) {
  th1 <- theta1_deg * pi / 180
  th2 <- theta2_deg * pi / 180
  th_outer <- seq(th1, th2, length.out = n)
  th_inner <- seq(th2, th1, length.out = n)

  tibble::tibble(
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

  if (nrow(plot_tbl) == 0) return(tibble::tibble())

  base_radius <- if (track == "outer") OUTER_TRACK_BASE else INNER_TRACK_BASE
  max_len     <- if (track == "outer") BAR_MAX_LEN_OUTER else BAR_MAX_LEN_INNER
  comp_col    <- if (track == "outer") "log2FC_A_vs_B" else "log2FC_C_vs_B"

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
        metabolite = rr$metabolite,
        region = rr$venn_region,
        comparison = ifelse(track == "outer", "GR vs MR", "MR+Phage cocktail vs MR"),
        direction = rr$direction,
        value = rr$comp_value,
        abs_value = rr$comp_abs,
        color = rr$color,
        poly_id = paste(track, rr$metabolite, sep = "__")
      )

    poly_list[[idx]] <- poly_df
    idx <- idx + 1L
  }

  dplyr::bind_rows(poly_list)
}

wrap_name <- function(x, width = 26) {
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
# 2. read metabolomics raw table
# -----------------------------
msg("Reading metabolomics raw table...")
stop_if_missing_file(METABOLOMICS_FILE)
met_raw <- readxl::read_excel(METABOLOMICS_FILE, sheet = METABOLOMICS_SHEET)
met_raw <- as.data.frame(met_raw)

required_met_cols <- c("Compound Name (nmol/g)", samples_keep)
if (!all(required_met_cols %in% colnames(met_raw))) {
  stop("Metabolomics table is missing required columns in sheet '", METABOLOMICS_SHEET, "'.")
}

meta_cols_met <- intersect(
  c("中文名称", "CAS", "Formula", "Pubchem.ID", "HMDB", "KEGG", "Class", "分类", "Average Molecular Weight"),
  colnames(met_raw)
)

met_long <- met_raw %>%
  select(all_of(c("Compound Name (nmol/g)", meta_cols_met, samples_keep))) %>%
  pivot_longer(cols = all_of(samples_keep), names_to = "sample", values_to = "value_raw") %>%
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

safe_write_csv(met_log2_by_sample, file.path(OUTDIR, "02_differential", "metabolite_log2_by_sample.csv"))

# -----------------------------
# 3. differential metabolites
# -----------------------------
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

    tibble::tibble(
      p_overall = p_overall,
      mean_A = mean_A, sd_A = sd_A, n_A = nA,
      mean_B = mean_B, sd_B = sd_B, n_B = nB,
      mean_C = mean_C, sd_C = sd_C, n_C = nC,
      log2FC_A_vs_B = log2FC_A_vs_B, p_A_vs_B = p_A_vs_B,
      log2FC_C_vs_B = log2FC_C_vs_B, p_C_vs_B = p_C_vs_B,
      log2FC_C_vs_A = log2FC_C_vs_A, p_C_vs_A = p_C_vs_A
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
    fdr_A_vs_B  = p.adjust(p_A_vs_B, method = "BH"),
    fdr_C_vs_B  = p.adjust(p_C_vs_B, method = "BH"),
    fdr_C_vs_A  = p.adjust(p_C_vs_A, method = "BH")
  ) %>%
  mutate(
    AB_sig = !is.na(fdr_A_vs_B) & fdr_A_vs_B < METAB_AB_FDR_CUT,
    CB_sig = !is.na(fdr_C_vs_B) & fdr_C_vs_B < METAB_CB_FDR_CUT,
    venn_region = case_when(
      AB_sig & CB_sig ~ "Shared",
      AB_sig & !CB_sig ~ "GR_vs_MR_only",
      !AB_sig & CB_sig ~ "Phage_vs_MR_only",
      TRUE ~ "Not_selected"
    ),
    order_key = pmax(abs(log2FC_A_vs_B), abs(log2FC_C_vs_B), na.rm = TRUE)
  ) %>%
  arrange(venn_region, desc(order_key), metabolite)

safe_write_csv(met_diff, file.path(OUTDIR, "02_differential", "metabolite_differential_all.csv"))

# -----------------------------
# 4. Venn sets and color map
# -----------------------------
msg("Building Venn sets...")

venn_tbl <- met_diff %>%
  filter(venn_region != "Not_selected") %>%
  mutate(region_order = case_when(
    venn_region == "GR_vs_MR_only" ~ 1L,
    venn_region == "Shared" ~ 2L,
    venn_region == "Phage_vs_MR_only" ~ 3L,
    TRUE ~ 99L
  )) %>%
  arrange(region_order, desc(order_key), metabolite)

# minimal modification:
# - only regions use venn-set colors
# - shared region metabolites use distinct colors
shared_metabolites <- venn_tbl %>%
  filter(venn_region == "Shared") %>%
  pull(metabolite) %>%
  unique()

shared_color_vec <- make_distinct_colors(length(shared_metabolites))
names(shared_color_vec) <- shared_metabolites

color_map_tbl <- venn_tbl %>%
  transmute(
    metabolite,
    venn_region,
    color = case_when(
      venn_region == "GR_vs_MR_only" ~ AONLY_COLOR,
      venn_region == "Phage_vs_MR_only" ~ BONLY_COLOR,
      venn_region == "Shared" ~ shared_color_vec[metabolite],
      TRUE ~ "grey80"
    )
  )

venn_tbl <- venn_tbl %>%
  left_join(color_map_tbl, by = c("metabolite", "venn_region"))

safe_write_csv(venn_tbl, file.path(OUTDIR, "03_venn_sets", "metabolite_venn_assignment_all.csv"))
safe_write_csv(color_map_tbl, file.path(OUTDIR, "03_venn_sets", "metabolite_color_map.csv"))
safe_write_csv(
  venn_tbl %>% filter(venn_region == "GR_vs_MR_only"),
  file.path(OUTDIR, "03_venn_sets", "GR_vs_MR_only.csv")
)
safe_write_csv(
  venn_tbl %>% filter(venn_region == "Shared"),
  file.path(OUTDIR, "03_venn_sets", "shared.csv")
)
safe_write_csv(
  venn_tbl %>% filter(venn_region == "Phage_vs_MR_only"),
  file.path(OUTDIR, "03_venn_sets", "Phage_vs_MR_only.csv")
)

# -----------------------------
# 5. assign angular layout
# -----------------------------
msg("Assigning ring layout...")

venn_A <- venn_tbl %>%
  filter(venn_region == "GR_vs_MR_only") %>%
  arrange(desc(order_key), metabolite)

venn_S <- venn_tbl %>%
  filter(venn_region == "Shared") %>%
  arrange(desc(order_key), metabolite)

venn_B <- venn_tbl %>%
  filter(venn_region == "Phage_vs_MR_only") %>%
  arrange(desc(order_key), metabolite)

venn_A <- assign_angles(venn_A, ANGLE_RANGE_AONLY[1], ANGLE_RANGE_AONLY[2])
venn_S <- assign_angles(venn_S, ANGLE_RANGE_SHARED[1], ANGLE_RANGE_SHARED[2])
venn_B <- assign_angles(venn_B, ANGLE_RANGE_BONLY[1], ANGLE_RANGE_BONLY[2])

plot_layout_tbl <- bind_rows(venn_A, venn_S, venn_B) %>%
  mutate(
    label = metabolite,
    label_wrapped = wrap_name(label, width = 24)
  )

safe_write_csv(plot_layout_tbl, file.path(OUTDIR, "04_plot", "ring_bar_layout.csv"))

# -----------------------------
# 6. build ring polygons
# -----------------------------
msg("Building ring bar polygons...")

max_abs_fc <- max(
  abs(c(plot_layout_tbl$log2FC_A_vs_B, plot_layout_tbl$log2FC_C_vs_B)),
  na.rm = TRUE
)
if (!is.finite(max_abs_fc) || max_abs_fc <= 0) max_abs_fc <- 1

outer_poly <- build_ring_poly(plot_layout_tbl, track = "outer", max_abs_fc = max_abs_fc)
inner_poly <- build_ring_poly(plot_layout_tbl, track = "inner", max_abs_fc = max_abs_fc)

ring_poly_tbl <- bind_rows(outer_poly, inner_poly)

safe_write_csv(ring_poly_tbl, file.path(OUTDIR, "04_plot", "ring_bar_segments.csv"))

# baseline circles for tracks
ring_base_outer <- make_circle_df(0, 0, OUTER_TRACK_BASE, n = 400)
ring_base_inner <- make_circle_df(0, 0, INNER_TRACK_BASE, n = 400)

venn_circle_left  <- make_circle_df(VENN_LEFT_X, VENN_CENTER_Y, VENN_RADIUS, n = 400)
venn_circle_right <- make_circle_df(VENN_RIGHT_X, VENN_CENTER_Y, VENN_RADIUS, n = 400)

# -----------------------------
# 7. plot
# -----------------------------
msg("Plotting Venn + outer ring bars...")

count_Aonly <- sum(plot_layout_tbl$venn_region == "GR_vs_MR_only")
count_Shared <- sum(plot_layout_tbl$venn_region == "Shared")
count_Bonly <- sum(plot_layout_tbl$venn_region == "Phage_vs_MR_only")

fill_values <- setNames(plot_layout_tbl$color, plot_layout_tbl$metabolite)

# legend: only show all shared metabolites
legend_df <- plot_layout_tbl %>%
  filter(venn_region == "Shared") %>%
  distinct(metabolite, label, color) %>%
  arrange(metabolite)

legend_breaks <- legend_df$metabolite
legend_labels <- legend_df$label

subtitle_text <- paste0(
  "Outer track = GR vs MR; inner track = MR+Phage cocktail vs MR. ",
  "Outward = higher in the first group of each comparison; inward = lower. ",
  "Bar length = |log2FC|."
)

p <- ggplot() +
  # venn circles
  geom_path(data = venn_circle_left, aes(x, y),
            color = "#7CAE00", linewidth = 0.8) +
  geom_path(data = venn_circle_right, aes(x, y),
            color = "#00BFC4", linewidth = 0.8) +

  # translucent fills for circles
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

  # ring baselines
  geom_path(data = ring_base_outer, aes(x, y), color = "grey85", linewidth = 0.4) +
  geom_path(data = ring_base_inner, aes(x, y), color = "grey88", linewidth = 0.4) +

  # ring polygons
  geom_polygon(
    data = ring_poly_tbl,
    aes(x = x, y = y, group = poly_id, fill = metabolite),
    color = NA,
    alpha = 0.95
  ) +

  # venn labels
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
    guide = guide_legend(ncol = 1, byrow = TRUE)
  ) +
  coord_equal(xlim = c(-2.35, 2.95), ylim = c(-2.20, 2.20), clip = "off") +
  labs(
    title = "Differential metabolite Venn with outer ring bars",
    subtitle = subtitle_text,
    fill = "Shared metabolites"
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 9),
    plot.margin = margin(12, 16, 12, 16)
  )

ggsave(
  filename = file.path(OUTDIR, "04_plot", "metabolite_venn_outer_ringbars.pdf"),
  plot = p, width = 14, height = 10
)
ggsave(
  filename = file.path(OUTDIR, "04_plot", "metabolite_venn_outer_ringbars.png"),
  plot = p, width = 14, height = 10, dpi = 300
)

msg("Done.")
msg("Output directory: ", normalizePath(OUTDIR, winslash = "/"))