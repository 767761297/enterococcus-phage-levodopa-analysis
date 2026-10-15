## ============================================================
## Figure 2I: Dopamine, final 3-group version
##
## Three-group plot for Dopamine from the metabolomics quantification table.
##
## Updated plotting rule:
##   - NO significance brackets/connecting lines
##   - Only place significance labels above:
##       GR
##       MR+Phage cocktail
##     to indicate comparison vs MR
##
## Significance labels:
##   ns    : p >= 0.05
##   *     : p < 0.05
##   **    : p < 0.01
##   ***   : p < 0.001
##   ****  : p < 0.0001
##
## Inputs:
##   - Q2.xlsx
##
## Outputs:
##   - per-sample raw data
##   - matched metabolite row info
##   - overall Kruskal-Wallis result
##   - pairwise Wilcoxon results
##   - violin plot PDF/PNG
## ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(readr)
  library(ggplot2)
})

# -------------------- SETTINGS --------------------
MET_QUANT_FILE <- "Q2.xlsx"
MET_SHEET <- NULL

# exact first, then fallback fuzzy "dopamin"
TARGET_METABOLITE_CANDIDATES <- c("Dopamine", "Dopamin")

OUTDIR  <- "results_Fig2I_Dopamine_3groups"
TAB_DIR <- file.path(OUTDIR, "tables")
PDF_DIR <- file.path(OUTDIR, "figures_pdf")
PNG_DIR <- file.path(OUTDIR, "figures_png")

for (d in c(OUTDIR, TAB_DIR, PDF_DIR, PNG_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

GROUP_LEVELS <- c("GR", "MR", "MR+Phage cocktail")
GROUP_COLORS <- c(
  "GR" = "#7CAE00",
  "MR" = "#F8766D",
  "MR+Phage cocktail" = "#00BFC4"
)

POINT_SIZE <- 2.5
POINT_SHAPE <- 16
VIOLIN_ALPHA <- 0.12
VIOLIN_WIDTH <- 0.78
BOX_WIDTH <- 0.18
JITTER_WIDTH <- 0.07

# -------------------- HELPERS --------------------
stopifnot(file.exists(MET_QUANT_FILE))
if (is.null(MET_SHEET)) MET_SHEET <- excel_sheets(MET_QUANT_FILE)[1]

map_group_3 <- function(s) {
  if (str_detect(s, "^A")) return("GR")
  if (str_detect(s, "^B")) return("MR")
  if (str_detect(s, "^C")) return("MR+Phage cocktail")
  return(NA_character_)
}

p_to_label <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-4) return("****")
  if (p < 1e-3) return("***")
  if (p < 1e-2) return("**")
  if (p < 5e-2) return("*")
  return("ns")
}

safe_wilcox <- function(x, y) {
  out <- tryCatch(
    suppressWarnings(wilcox.test(x, y, exact = FALSE)),
    error = function(e) NULL
  )
  if (is.null(out)) {
    return(list(statistic = NA_real_, p.value = NA_real_))
  }
  list(statistic = unname(out$statistic), p.value = out$p.value)
}

theme_violin_clean <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      legend.position = "none"
    )
}

# -------------------- READ INPUT --------------------
met_raw <- read_excel(MET_QUANT_FILE, sheet = MET_SHEET) %>%
  as.data.frame(check.names = FALSE)

sample_cols <- grep("^[ABC][1-5]$", names(met_raw), value = TRUE)
if (length(sample_cols) != 15) {
  stop("Cannot detect 15 sample columns A1..C5 in Q2.xlsx. Found: ",
       paste(sample_cols, collapse = ", "))
}

compound_col <- grep("Compound\\s*Name.*nmol/g", names(met_raw), ignore.case = TRUE, value = TRUE)[1]
if (is.na(compound_col) || compound_col == "") {
  compound_col <- grep("Compound\\s*Name", names(met_raw), ignore.case = TRUE, value = TRUE)[1]
}
if (is.na(compound_col) || compound_col == "") {
  stop("Cannot find 'Compound Name (nmol/g)' or 'Compound Name' column.")
}

met_name <- met_raw[[compound_col]] %>%
  as.character() %>%
  str_trim()

# -------------------- MATCH TARGET METABOLITE --------------------
target_hit_exact <- which(tolower(met_name) %in% tolower(TARGET_METABOLITE_CANDIDATES))
if (length(target_hit_exact) == 0) {
  target_hit_exact <- grep("dopamin", met_name, ignore.case = TRUE)
}
if (length(target_hit_exact) == 0) {
  stop("Cannot find Dopamine/Dopamin in Q2.xlsx. Please check the exact metabolite name.")
}

matched_rows <- met_raw[target_hit_exact, , drop = FALSE]
matched_names <- met_name[target_hit_exact]

write_csv(
  tibble(row_index = target_hit_exact, matched_metabolite_name = matched_names),
  file.path(TAB_DIR, "matched_metabolite_rows.csv")
)

# if duplicate matched rows exist, average by sample
dopamine_mat <- matched_rows[, sample_cols, drop = FALSE] %>% as.matrix()
dopamine_mat <- apply(dopamine_mat, 2, function(x) suppressWarnings(as.numeric(x)))

if (is.null(dim(dopamine_mat))) {
  dopamine_values <- as.numeric(dopamine_mat)
  names(dopamine_values) <- sample_cols
} else {
  dopamine_values <- colMeans(dopamine_mat, na.rm = TRUE)
  names(dopamine_values) <- sample_cols
}

# -------------------- BUILD PER-SAMPLE TABLE --------------------
df <- tibble(
  sample = sample_cols,
  group = factor(vapply(sample_cols, map_group_3, character(1)), levels = GROUP_LEVELS),
  Dopamine = as.numeric(dopamine_values[sample_cols])
) %>%
  arrange(factor(sample, levels = sample_cols))

write_csv(df, file.path(TAB_DIR, "Dopamine_per_sample_raw.csv"))
write_tsv(df, file.path(TAB_DIR, "Dopamine_per_sample_raw.tsv"))

# -------------------- STATS --------------------
kw_p <- kruskal.test(Dopamine ~ group, data = df)$p.value
kw_tbl <- tibble(
  metabolite = "Dopamine",
  test = "Kruskal-Wallis",
  p_value = kw_p
)
write_csv(kw_tbl, file.path(TAB_DIR, "Dopamine_overall_KruskalWallis.csv"))

pairs <- combn(GROUP_LEVELS, 2, simplify = FALSE)
pw_tbl <- lapply(pairs, function(gp) {
  g1 <- gp[1]
  g2 <- gp[2]
  x1 <- df %>% filter(group == g1) %>% pull(Dopamine)
  x2 <- df %>% filter(group == g2) %>% pull(Dopamine)
  wt <- safe_wilcox(x1, x2)

  tibble(
    metabolite = "Dopamine",
    group1 = g1,
    group2 = g2,
    n1 = length(x1),
    n2 = length(x2),
    mean1 = mean(x1, na.rm = TRUE),
    mean2 = mean(x2, na.rm = TRUE),
    median1 = median(x1, na.rm = TRUE),
    median2 = median(x2, na.rm = TRUE),
    delta_mean = mean1 - mean2,
    delta_median = median1 - median2,
    statistic_W = as.numeric(wt$statistic),
    p_value_raw = as.numeric(wt$p.value),
    signif_label_raw = p_to_label(wt$p.value)
  )
}) %>% bind_rows()

pw_tbl$p_value_BH <- p.adjust(pw_tbl$p_value_raw, method = "BH")

write_csv(pw_tbl, file.path(TAB_DIR, "Dopamine_pairwise_Wilcoxon_rawP_BH.csv"))
write_tsv(pw_tbl, file.path(TAB_DIR, "Dopamine_pairwise_Wilcoxon_rawP_BH.tsv"))

# -------------------- PLOT LABELS --------------------
# Final Figure 2I includes three groups only. Display the two comparisons versus MR.
label_tbl <- bind_rows(
  pw_tbl %>% filter((group1 == "GR" & group2 == "MR") | (group1 == "MR" & group2 == "GR")) %>%
    slice(1) %>% mutate(target_group = "GR"),
  pw_tbl %>% filter((group1 == "MR+Phage cocktail" & group2 == "MR") |
                    (group1 == "MR" & group2 == "MR+Phage cocktail")) %>%
    slice(1) %>% mutate(target_group = "MR+Phage cocktail")
)

y_max <- max(df$Dopamine, na.rm = TRUE)
y_min <- min(df$Dopamine, na.rm = TRUE)
y_range <- y_max - y_min
if (!is.finite(y_range) || y_range <= 0) y_range <- max(abs(y_max), 1)

offset1 <- max(0.08 * y_range, 0.8)
offset2 <- max(0.16 * y_range, 1.6)

label_pos <- tibble(
  target_group = c("GR", "MR+Phage cocktail"),
  x = c(1, 3),
  y = c(
    max(df$Dopamine[df$group == "GR"], na.rm = TRUE) + offset1,
    max(df$Dopamine[df$group == "MR+Phage cocktail"], na.rm = TRUE) + offset2
  )
)

label_tbl <- label_tbl %>%
  left_join(label_pos, by = "target_group")

# -------------------- PLOT --------------------
p <- ggplot(df, aes(x = group, y = Dopamine, fill = group, color = group)) +
  geom_violin(alpha = VIOLIN_ALPHA, trim = TRUE, scale = "width", width = VIOLIN_WIDTH, linewidth = 0.8) +
  geom_boxplot(width = BOX_WIDTH, outlier.shape = NA, fill = NA, linewidth = 0.8) +
  geom_point(
    position = position_jitter(width = JITTER_WIDTH, height = 0),
    size = POINT_SIZE,
    shape = POINT_SHAPE,
    stroke = 0
  ) +
  geom_text(
    data = label_tbl,
    aes(x = x, y = y, label = signif_label_raw),
    inherit.aes = FALSE,
    size = 6,
    fontface = "bold"
  ) +
  scale_color_manual(values = GROUP_COLORS, drop = FALSE) +
  scale_fill_manual(values = GROUP_COLORS, drop = FALSE) +
  scale_x_discrete(labels = c("GR", "MR", "MR + Phage")) +
  labs(
    title = "Dopamine",
    x = NULL,
    y = "Dopamine (nmol/g)"
  ) +
  theme_violin_clean() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

fig_w <- 6.2
fig_h <- 6.0

ggsave(file.path(PDF_DIR, "Dopamine_3groups_vsMR_labelOnly.pdf"), p, width = fig_w, height = fig_h)
ggsave(file.path(PNG_DIR, "Dopamine_3groups_vsMR_labelOnly.png"), p, width = fig_w, height = fig_h, dpi = 300)

message("[DONE] Outputs saved in: ", normalizePath(OUTDIR))