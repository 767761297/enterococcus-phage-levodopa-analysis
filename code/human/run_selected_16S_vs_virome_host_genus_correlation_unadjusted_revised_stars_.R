## =========================================
## Selected 16S genus vs selected virome host-genus correlation bubble plot
## Unadjusted / CLR / Spearman / BH / revised stars
##
## Required input files in working directory:
##   metadata.csv
##   table_f_even_depth_w_tax.xls
##   vOTU_abundance_matrix_TPM.tsv
##   cherry_prediction.tsv
##
## Star rule:
##   *   : p < 0.05 and q > 0.1
##   **  : p < 0.05 and q < 0.1
##   *** : p < 0.05 and q < 0.05
##   q = FDR (BH-adjusted q value)
## =========================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(tibble)
})

## ---------- inputs ----------
META_FILE  <- "metadata.csv"
INPUT_16S  <- "table_filtered_w_tax.xls"
TPM_FILE   <- "vOTU_abundance_matrix_TPM.tsv"
HOST_FILE  <- "cherry_prediction.tsv"

OUT_ROOT <- "results_16S_virome_bubbleplot"
dir.create(OUT_ROOT, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "00_manifest"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "01_feature_lists"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "02_per_sample_matrices"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "03_stats"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "04_figures_pdf"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "04_figures_png"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "05_scripts"), recursive = TRUE, showWarnings = FALSE)

## ---------- params ----------
GROUP_LEVELS <- c("PD_GR", "PD_MR")
PSEUDOCOUNT_16S <- 0.5
PSEUDOCOUNT_VIROME <- 0.01
MIN_PREVALENCE_16S <- 0.10
MIN_TOTAL_COUNT_16S <- 10
MIN_PREVALENCE_VOTU <- 0.10

FEATURES_16S <- c(
  "Enterococcus","Solobacterium","Carnobacterium","Veillonella","Streptococcus",
  "Xylanibacter","Serratia","Olivibacter","Niallia","Dysosmobacter",
  "Acidaminococcus","Achromobacter","Anaerotruncus"
)

FEATURES_HOST_GENUS <- c(
  "Bacillus","Sphingobium","[Ruminococcus]","Herbinix","Bittarella",
  "Lactonifactor","Burkholderiales","Oceanobacillus","Leptospira",
  "Enterococcus","Raoultella","Salinispora","Bifidobacterium",
  "Odoribacter","Bacteroidetes"
)

FEATURES_HOST_DISPLAY <- c(
  "Bacillus | top:Bacillus subtilis",
  "Sphingobium | top:Sphingobium sp. IP1",
  "[Ruminococcus] | top:[Ruminococcus] gnavus",
  "Herbinix | top:Herbinix hemicellulosilytica",
  "Bittarella | top:Bittarella massiliensis",
  "Lactonifactor | top:Lactonifactor longoviformis",
  "Burkholderiales | top:Burkholderiales bacterium",
  "Oceanobacillus | top:Oceanobacillus caeni",
  "Leptospira | top:Leptospira mayottensis",
  "Enterococcus | top:Enterococcus faecalis",
  "Raoultella | top:Raoultella planticola",
  "Salinispora | top:Salinispora arenicola",
  "Bifidobacterium | top:Bifidobacterium moukalabense",
  "Odoribacter | top:Odoribacter splanchnicus",
  "Bacteroidetes | top:Bacteroidetes bacterium oral taxon 272"
)

## ---------- helpers ----------
read_csv_smart <- function(path) {
  tryCatch(
    readr::read_csv(path, show_col_types = FALSE, locale = locale(encoding = "GB18030")),
    error = function(e) readr::read_csv(path, show_col_types = FALSE, locale = locale(encoding = "UTF-8"))
  )
}

normalize_id <- function(x) {
  x |>
    as.character() |>
    str_trim() |>
    str_replace_all("\\s+", "")
}

write_tsv0 <- function(x, path) readr::write_tsv(x, path)

parse_tax_16s <- function(tax) {
  parts <- str_split(as.character(tax), ";", simplify = TRUE)
  get_rank <- function(row, key) {
    cand <- row[str_detect(row, paste0("^\\s*", key, "_+"))]
    if (length(cand) == 0) return(NA_character_)
    v <- str_trim(cand[1])
    v <- str_replace(v, paste0("^\\s*", key, "_+"), "")
    if (v == "" || v == "NA") NA_character_ else v
  }
  tibble(genus = apply(parts, 1, get_rank, key = "g"))
}

fill_unknown <- function(x, label) {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "NA"] <- paste0("Unassigned_", label)
  x
}

aggregate_counts <- function(mat_feat_by_sample, feature_to_group_named) {
  grp <- feature_to_group_named[rownames(mat_feat_by_sample)]
  rowsum(mat_feat_by_sample, group = grp, reorder = TRUE)
}

clr_transform <- function(mat_feat_by_sample, pseudocount) {
  x <- mat_feat_by_sample + pseudocount
  x <- sweep(x, 2, colSums(x), "/")
  lx <- log(x)
  sweep(lx, 2, colMeans(lx), "-")
}

parse_host_display <- function(host_field) {
  host_field <- as.character(host_field)
  host_field[host_field %in% c("-", "", "NA")] <- NA_character_

  rank <- ifelse(!is.na(host_field) & str_detect(host_field, ":"), str_split_fixed(host_field, ":", 2)[, 1], NA_character_)
  name <- ifelse(!is.na(host_field) & str_detect(host_field, ":"), str_split_fixed(host_field, ":", 2)[, 2], host_field)
  name <- str_trim(name)

  host_species <- ifelse(!is.na(rank) & rank == "species" & !is.na(name) & name != "", name, NA_character_)
  host_genus <- ifelse(!is.na(rank) & rank == "genus" & !is.na(name) & name != "", name,
                       ifelse(!is.na(host_species) & str_detect(host_species, "\\s+"),
                              str_split_fixed(host_species, "\\s+", 2)[, 1],
                              NA_character_))
  host_genus <- ifelse(!is.na(host_genus) & host_genus != "", host_genus, "Unassigned_host")

  tibble(host_genus = host_genus, host_species = host_species)
}

safe_mode <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != "" & x != "NA"]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

stop_if_missing_cols <- function(df, cols, file_label) {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) {
    stop(sprintf("%s missing required column(s): %s", file_label, paste(miss, collapse = ", ")))
  }
}

## ---------- metadata ----------
meta_raw <- read_csv_smart(META_FILE)
stop_if_missing_cols(meta_raw, c("Sample", "Group"), META_FILE)

meta <- meta_raw |>
  transmute(
    sample_id = normalize_id(Sample),
    group = Group
  ) |>
  filter(group %in% GROUP_LEVELS) |>
  mutate(group = factor(group, levels = GROUP_LEVELS)) |>
  arrange(group, sample_id)

sample_order <- meta$sample_id
if (length(sample_order) == 0) stop("No PD_GR / PD_MR samples found in metadata.csv")

## ---------- 16S ----------
dt16 <- fread(INPUT_16S, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
if (!"taxonomy" %in% names(dt16)) stop("table_f_even_depth_w_tax.xls must contain a 'taxonomy' column")

sample_cols_16s_missing <- setdiff(sample_order, names(dt16))
if (length(sample_cols_16s_missing) > 0) {
  stop(sprintf("16S table missing sample columns: %s", paste(sample_cols_16s_missing, collapse = ", ")))
}

asv_id <- as.character(dt16[[1]])
mat_asv <- dt16[, sample_order, drop = FALSE]
mat_asv[] <- lapply(mat_asv, as.numeric)
mat_asv <- as.matrix(mat_asv)
rownames(mat_asv) <- asv_id

tax_tbl <- tibble(asv_id = asv_id, taxonomy_raw = dt16$taxonomy) |>
  bind_cols(parse_tax_16s(dt16$taxonomy)) |>
  mutate(genus = fill_unknown(genus, "genus"))

map_asv_genus <- setNames(tax_tbl$genus, tax_tbl$asv_id)
mat_genus <- aggregate_counts(mat_asv, map_asv_genus)
prev_all <- rowMeans(mat_genus > 0)
tot <- rowSums(mat_genus)
keep16 <- (prev_all >= MIN_PREVALENCE_16S) & (tot >= MIN_TOTAL_COUNT_16S)
mat_genus_f <- mat_genus[keep16, , drop = FALSE]

missing_features_16s <- setdiff(FEATURES_16S, rownames(mat_genus_f))
if (length(missing_features_16s) > 0) {
  stop(sprintf("Selected 16S genus feature(s) absent after filtering: %s", paste(missing_features_16s, collapse = ", ")))
}

mat16_raw <- mat_genus_f[FEATURES_16S, sample_order, drop = FALSE]
mat16_clr <- clr_transform(mat16_raw, PSEUDOCOUNT_16S)

## ---------- virome host-genus ----------
tpm_dt <- fread(TPM_FILE, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
if (!"id" %in% names(tpm_dt)) stop("vOTU_abundance_matrix_TPM.tsv must contain an 'id' column")

sample_cols_v_missing <- setdiff(sample_order, names(tpm_dt))
if (length(sample_cols_v_missing) > 0) {
  stop(sprintf("Virome TPM table missing sample columns: %s", paste(sample_cols_v_missing, collapse = ", ")))
}

mat_votu <- tpm_dt[, sample_order, drop = FALSE]
mat_votu[] <- lapply(mat_votu, as.numeric)
mat_votu <- as.matrix(mat_votu)
rownames(mat_votu) <- as.character(tpm_dt$id)

prev_votu <- rowMeans(mat_votu > 0)
keep_votu <- prev_votu >= MIN_PREVALENCE_VOTU
mat_votu0 <- mat_votu[keep_votu, , drop = FALSE]

host_dt <- fread(HOST_FILE, sep = "\t", header = TRUE, data.table = FALSE)
stop_if_missing_cols(host_dt, c("Accession", "Host"), HOST_FILE)

host_tbl <- host_dt |>
  transmute(
    votu_id = as.character(Accession),
    Host_raw = as.character(Host)
  ) |>
  bind_cols(parse_host_display(host_dt$Host))

map_host <- setNames(host_tbl$host_genus, host_tbl$votu_id)
map_host <- map_host[intersect(names(map_host), rownames(mat_votu0))]
if (length(map_host) == 0) stop("No overlapping vOTU IDs between virome TPM table and cherry_prediction.tsv")

mat_host <- aggregate_counts(mat_votu0[names(map_host), , drop = FALSE], map_host)
missing_features_host <- setdiff(FEATURES_HOST_GENUS, rownames(mat_host))
if (length(missing_features_host) > 0) {
  stop(sprintf("Selected virome host-genus feature(s) absent after aggregation/filtering: %s", paste(missing_features_host, collapse = ", ")))
}

mat_host_raw <- mat_host[FEATURES_HOST_GENUS, sample_order, drop = FALSE]
mat_host_clr <- clr_transform(mat_host_raw, PSEUDOCOUNT_VIROME)

display_tbl <- host_tbl |>
  filter(votu_id %in% rownames(mat_votu0)) |>
  filter(host_genus %in% FEATURES_HOST_GENUS) |>
  group_by(host_genus) |>
  summarise(top_host_species = safe_mode(host_species), .groups = "drop") |>
  mutate(display_label = ifelse(!is.na(top_host_species) & top_host_species != "",
                                paste0(host_genus, " | top:", top_host_species),
                                host_genus))

## Force full-label order consistent with final heatmap
host_display_map <- setNames(FEATURES_HOST_DISPLAY, FEATURES_HOST_GENUS)
display_tbl <- tibble(
  host_genus = FEATURES_HOST_GENUS,
  display_label = FEATURES_HOST_DISPLAY
) |>
  left_join(display_tbl |> select(host_genus, top_host_species), by = "host_genus")

## ---------- exports ----------
write_tsv0(meta, file.path(OUT_ROOT, "00_manifest", "metadata_40_samples.tsv"))
write_tsv0(tibble(feature_16S_genus = FEATURES_16S), file.path(OUT_ROOT, "01_feature_lists", "selected_16S_genus_features.tsv"))
write_tsv0(display_tbl, file.path(OUT_ROOT, "01_feature_lists", "selected_virome_host_genus_features.tsv"))

write_tsv0(as.data.frame(t(mat16_raw)) |> rownames_to_column("sample_id"),
           file.path(OUT_ROOT, "02_per_sample_matrices", "selected_16S_genus_raw_counts_samples_x_features.tsv"))
write_tsv0(as.data.frame(t(mat16_clr)) |> rownames_to_column("sample_id"),
           file.path(OUT_ROOT, "02_per_sample_matrices", "selected_16S_genus_CLR_samples_x_features.tsv"))

host_raw_export <- as.data.frame(t(mat_host_raw))
colnames(host_raw_export) <- FEATURES_HOST_DISPLAY
host_clr_export <- as.data.frame(t(mat_host_clr))
colnames(host_clr_export) <- FEATURES_HOST_DISPLAY
write_tsv0(host_raw_export |> rownames_to_column("sample_id"),
           file.path(OUT_ROOT, "02_per_sample_matrices", "selected_virome_host_genus_raw_TPM_samples_x_features.tsv"))
write_tsv0(host_clr_export |> rownames_to_column("sample_id"),
           file.path(OUT_ROOT, "02_per_sample_matrices", "selected_virome_host_genus_CLR_samples_x_features.tsv"))

## ---------- correlations ----------
res <- list()
for (g16 in FEATURES_16S) {
  x <- as.numeric(mat16_clr[g16, ])
  for (hg in FEATURES_HOST_GENUS) {
    y <- as.numeric(mat_host_clr[hg, ])
    tmp <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
    res[[length(res) + 1]] <- tibble(
      feature_16S_genus = g16,
      feature_virome_host_genus = hg,
      feature_virome_display = host_display_map[[hg]],
      rho = unname(tmp$estimate),
      p_value = tmp$p.value,
      n_samples = length(x)
    )
  }
}

res_df <- bind_rows(res) |>
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
    star = case_when(
      p_value < 0.05 & q_value < 0.05 ~ "***",
      p_value < 0.05 & q_value < 0.10 ~ "**",
      p_value < 0.05 & q_value > 0.10 ~ "*",
      TRUE ~ ""
    ),
    abs_rho = abs(rho),
    feature_16S_genus = factor(feature_16S_genus, levels = FEATURES_16S),
    feature_virome_display = factor(feature_virome_display, levels = FEATURES_HOST_DISPLAY),
    feature_virome_host_genus = factor(feature_virome_host_genus, levels = FEATURES_HOST_GENUS)
  )

write_tsv0(res_df, file.path(OUT_ROOT, "03_stats", "pairwise_spearman_CLR_16S_vs_virome_host_genus_long_revised_stars.tsv"))

rho_mat <- res_df |>
  select(feature_16S_genus, feature_virome_display, rho) |>
  pivot_wider(names_from = feature_virome_display, values_from = rho)
p_mat <- res_df |>
  select(feature_16S_genus, feature_virome_display, p_value) |>
  pivot_wider(names_from = feature_virome_display, values_from = p_value)
q_mat <- res_df |>
  select(feature_16S_genus, feature_virome_display, q_value) |>
  pivot_wider(names_from = feature_virome_display, values_from = q_value)
star_mat_full <- res_df |>
  select(feature_16S_genus, feature_virome_display, star) |>
  pivot_wider(names_from = feature_virome_display, values_from = star)
star_mat_short <- res_df |>
  select(feature_16S_genus, feature_virome_host_genus, star) |>
  pivot_wider(names_from = feature_virome_host_genus, values_from = star)

write_tsv0(rho_mat, file.path(OUT_ROOT, "03_stats", "rho_matrix_spearman_CLR.tsv"))
write_tsv0(p_mat, file.path(OUT_ROOT, "03_stats", "p_matrix_spearman_CLR.tsv"))
write_tsv0(q_mat, file.path(OUT_ROOT, "03_stats", "q_matrix_spearman_CLR.tsv"))
write_tsv0(star_mat_full, file.path(OUT_ROOT, "03_stats", "star_matrix_revised_full_labels.tsv"))
write_tsv0(star_mat_short, file.path(OUT_ROOT, "03_stats", "star_matrix_revised_short_labels.tsv"))

## ---------- plotting data ----------
plot_df_short <- res_df |>
  mutate(
    x = factor(feature_virome_host_genus, levels = FEATURES_HOST_GENUS),
    y = factor(feature_16S_genus, levels = rev(FEATURES_16S))
  )

plot_df_full <- res_df |>
  mutate(
    x = factor(feature_virome_display, levels = FEATURES_HOST_DISPLAY),
    y = factor(feature_16S_genus, levels = rev(FEATURES_16S))
  )

common_theme <- theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_line(color = "grey90")
  )

subtitle_txt <- "* p < 0.05 and q > 0.1; ** p < 0.05 and q < 0.1; *** p < 0.05 and q < 0.05 (q = FDR)"

p_short <- ggplot(plot_df_short, aes(x = x, y = y)) +
  geom_point(aes(size = abs_rho, fill = rho), shape = 21, color = "grey50", stroke = 0.3) +
  geom_text(aes(label = star), size = 3, fontface = "bold") +
  scale_fill_gradient2(low = "#3B6BD6", mid = "white", high = "#D81B36", midpoint = 0, limits = c(-0.6, 0.6)) +
  scale_size_area(max_size = 14, breaks = c(0.2, 0.4, 0.6)) +
  common_theme +
  labs(
    x = "Virome host-genus features",
    y = "16S genus features",
    fill = "Spearman's rho",
    size = "|rho|",
    title = "Selected 16S genera vs selected virome host-genera (CLR, unadjusted)",
    subtitle = subtitle_txt
  )

p_full <- ggplot(plot_df_full, aes(x = x, y = y)) +
  geom_point(aes(size = abs_rho, fill = rho), shape = 21, color = "grey50", stroke = 0.3) +
  geom_text(aes(label = star), size = 3, fontface = "bold") +
  scale_fill_gradient2(low = "#3B6BD6", mid = "white", high = "#D81B36", midpoint = 0, limits = c(-0.6, 0.6)) +
  scale_size_area(max_size = 14, breaks = c(0.2, 0.4, 0.6)) +
  common_theme +
  labs(
    x = "Virome host-genus features",
    y = "16S genus features",
    fill = "Spearman's rho",
    size = "|rho|",
    title = "Selected 16S genera vs selected virome host-genera (CLR, unadjusted)",
    subtitle = subtitle_txt
  )

ggsave(file.path(OUT_ROOT, "04_figures_pdf", "bubbleplot_selected_16S_vs_virome_host_genus_spearman_CLR_unadjusted_revised_stars_shortlabels.pdf"),
       p_short, width = 12.5, height = 8.8)
ggsave(file.path(OUT_ROOT, "04_figures_png", "bubbleplot_selected_16S_vs_virome_host_genus_spearman_CLR_unadjusted_revised_stars_shortlabels.png"),
       p_short, width = 12.5, height = 8.8, dpi = 300)

ggsave(file.path(OUT_ROOT, "04_figures_pdf", "bubbleplot_selected_16S_vs_virome_host_genus_spearman_CLR_unadjusted_revised_stars.pdf"),
       p_full, width = 17.0, height = 8.8)
ggsave(file.path(OUT_ROOT, "04_figures_png", "bubbleplot_selected_16S_vs_virome_host_genus_spearman_CLR_unadjusted_revised_stars.png"),
       p_full, width = 17.0, height = 8.8, dpi = 300)

## ---------- save script copy ----------
script_path <- tryCatch(normalizePath(sys.frames()[[1]]$ofile), error = function(e) NA_character_)
if (!is.na(script_path) && file.exists(script_path)) {
  file.copy(script_path,
            file.path(OUT_ROOT, "05_scripts", basename(script_path)),
            overwrite = TRUE)
}

cat("========================================\n")
cat("Done.\n")
cat("Output directory:\n")
cat(normalizePath(OUT_ROOT), "\n")
cat("========================================\n")
