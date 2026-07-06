## =========================
## 02_beta_16S_ASV_full.R  (FULL FIX: group join everywhere)
## =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(readr)
  library(janitor)
  library(vegan)
  library(ggplot2)
  library(ape)
  library(phyloseq)
  library(tibble)
})

GROUP_LEVELS <- c("PD_GR", "PD_MR")
GROUP_COLORS <- c(PD_GR = "#7CAE00", PD_MR = "#F8766D")

## ---------- helpers ----------
detect_file <- function(prefer, patterns) {
  if (!is.null(prefer) && file.exists(prefer)) return(prefer)
  files <- list.files(".", full.names = FALSE)
  for (pat in patterns) {
    hit <- files[str_detect(files, regex(pat, ignore_case = TRUE))]
    if (length(hit) > 0) return(hit[1])
  }
  stop("Cannot find input file. Patterns: ", paste(patterns, collapse = " | "))
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

dir_init <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(path, "tables"), showWarnings = FALSE)
  dir.create(file.path(path, "figures"), showWarnings = FALSE)
  dir.create(file.path(path, "logs"), showWarnings = FALSE)
  dir.create(file.path(path, "00_manifest"), showWarnings = FALSE)
}

write_tsv0 <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(x, path)
}

plot_scatter_ellipse <- function(df, x, y, title, xlab, ylab, out_pdf, out_png) {
  df <- df |> mutate(group = factor(group, levels = GROUP_LEVELS))

  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]], color = group, fill = group)) +
    stat_ellipse(aes(group = group, fill = group),
                 geom = "polygon", type = "t", level = 0.95,
                 alpha = 0.20, color = NA) +
    stat_ellipse(aes(group = group, color = group),
                 geom = "path", type = "t", level = 0.95,
                 linewidth = 0.8) +
    geom_point(size = 3) +
    scale_color_manual(values = GROUP_COLORS, breaks = GROUP_LEVELS, drop = FALSE) +
    scale_fill_manual(values = GROUP_COLORS, breaks = GROUP_LEVELS, drop = FALSE) +
    theme_bw() +
    labs(title = title, x = xlab, y = ylab, color = "group", fill = "group")

  pdf(out_pdf, width = 6.8, height = 5.2)
  suppressWarnings(print(p)); dev.off()

  png(out_png, width = 1800, height = 1400, res = 300)
  suppressWarnings(print(p)); dev.off()
}

write_dist_tables <- function(dist_obj, out_prefix) {
  mat <- as.matrix(dist_obj)
  write_tsv0(as.data.frame(mat) |> rownames_to_column("sample_id"),
             paste0(out_prefix, "_distance_matrix_wide.tsv"))
  long <- as.data.frame(as.table(mat), stringsAsFactors = FALSE) |>
    rename(sample1 = Var1, sample2 = Var2, distance = Freq) |>
    filter(sample1 < sample2)
  write_tsv0(long, paste0(out_prefix, "_distance_matrix_long.tsv"))
}

pcoa_with_correction <- function(dist_obj) {
  p0 <- ape::pcoa(dist_obj, correction = "none")
  neg <- sum(p0$values$Eigenvalues < -1e-8, na.rm = TRUE)
  if (neg > 0) {
    p1 <- ape::pcoa(dist_obj, correction = "cailliez")
    attr(p1, "pcoa_correction") <- "cailliez"
    attr(p1, "neg_eig_count_raw") <- neg
    p1
  } else {
    attr(p0, "pcoa_correction") <- "none"
    attr(p0, "neg_eig_count_raw") <- 0L
    p0
  }
}

extract_rel_eig <- function(pcoa_obj) {
  eig <- pcoa_obj$values$Eigenvalues
  rel <- pcoa_obj$values$Relative_eig
  if (is.null(rel) || length(rel) == 0) rel <- pcoa_obj$values$Rel_corr_eig
  if (is.null(rel) || length(rel) == 0) {
    pos <- eig; pos[pos < 0] <- 0
    rel <- pos / sum(pos)
  }
  n <- min(length(eig), length(rel))
  list(eig = eig[seq_len(n)], rel = rel[seq_len(n)])
}

make_meta_map <- function(meta_df) {
  meta_df |> as.data.frame() |>
    rownames_to_column("sample_id") |>
    mutate(sample_id = normalize_id(sample_id),
           group = factor(as.character(group), levels = GROUP_LEVELS))
}

attach_meta <- function(df, meta_map, out_dir = NULL, tag = NULL, ...) {
  out <- df |>
    mutate(sample_id = normalize_id(sample_id)) |>
    left_join(meta_map, by = "sample_id")

  if (!is.null(out_dir) && !is.null(tag)) {
    miss <- out |> filter(is.na(group)) |> select(sample_id) |> distinct()
    if (nrow(miss) > 0) {
      write_tsv0(miss, file.path(out_dir, "tables", paste0("join_missing_group_", tag, ".tsv")))
    }
  }
  out
}

run_adonis_suite <- function(dist_obj, meta_df, out_dir, perm = 999) {
  df <- meta_df |> as.data.frame()

  a0 <- adonis2(dist_obj ~ group, data = df, permutations = perm, by = "margin")
  write_tsv0(as.data.frame(a0) |> rownames_to_column("term"),
             file.path(out_dir, "tables", "adonis2_unadjusted.tsv"))

  a1 <- adonis2(dist_obj ~ group + age + bmi + disease_duration + constipation,
                data = df, permutations = perm, by = "margin")
  write_tsv0(as.data.frame(a1) |> rownames_to_column("term"),
             file.path(out_dir, "tables", "adonis2_adjusted_primary.tsv"))

  if (!all(is.na(df$levodopa_dose))) {
    a2 <- adonis2(dist_obj ~ group + age + bmi + disease_duration + constipation + levodopa_dose,
                  data = df, permutations = perm, by = "margin")
    write_tsv0(as.data.frame(a2) |> rownames_to_column("term"),
               file.path(out_dir, "tables", "adonis2_adjusted_sensitivity_add_dose.tsv"))
  }
}

run_betadisper_suite <- function(dist_obj, meta_df, out_dir, perm = 999) {
  meta_map <- make_meta_map(meta_df)

  bd <- betadisper(dist_obj, group = meta_df$group)

  per_raw <- tibble(
    sample_id = names(bd$distances),
    distance_to_centroid = as.numeric(bd$distances)
  )

  per <- attach_meta(per_raw, meta_map, out_dir = out_dir, tag = "betadisper")

  write_tsv0(per, file.path(out_dir, "tables", "betadisper_distance_to_centroid_per_sample.tsv"))

  bd_aov <- anova(bd)
  write_tsv0(as.data.frame(bd_aov) |> rownames_to_column("term"),
             file.path(out_dir, "tables", "betadisper_anova.tsv"))

  bd_perm <- permutest(bd, permutations = perm)
  write_tsv0(as.data.frame(bd_perm$tab) |> rownames_to_column("term"),
             file.path(out_dir, "tables", "betadisper_permutest.tsv"))

  p <- ggplot(per, aes(x = group, y = distance_to_centroid, color = group, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.2) +
    geom_jitter(width = 0.15, height = 0, size = 2) +
    scale_color_manual(values = GROUP_COLORS, breaks = GROUP_LEVELS, drop = FALSE) +
    scale_fill_manual(values = GROUP_COLORS, breaks = GROUP_LEVELS, drop = FALSE) +
    theme_bw() +
    labs(x = NULL, y = "Distance to centroid", color = "group", fill = "group")

  pdf(file.path(out_dir, "figures", "betadisper_distance_to_centroid_boxplot.pdf"), width = 6.5, height = 4.8)
  suppressWarnings(print(p)); dev.off()

  png(file.path(out_dir, "figures", "betadisper_distance_to_centroid_boxplot.png"), width = 1800, height = 1300, res = 300)
  suppressWarnings(print(p)); dev.off()
}

run_ordination_outputs <- function(dist_obj, meta_df, out_dir, title_prefix) {
  meta_map <- make_meta_map(meta_df)

  p <- pcoa_with_correction(dist_obj)
  corr <- attr(p, "pcoa_correction")
  neg  <- attr(p, "neg_eig_count_raw")

  coords_raw <- as.data.frame(p$vectors[, 1:6, drop = FALSE]) |>
    rownames_to_column("sample_id")
  coords <- attach_meta(coords_raw, meta_map, out_dir = out_dir, tag = "pcoa")

  ev <- extract_rel_eig(p)
  eig <- tibble(
    axis = paste0("PCo", seq_along(ev$eig)),
    eigenvalue = ev$eig,
    rel_eig = ev$rel,
    cum_rel = cumsum(ev$rel)
  )

  write_tsv0(coords, file.path(out_dir, "tables", "pcoa_scores_per_sample.tsv"))
  write_tsv0(eig,    file.path(out_dir, "tables", "pcoa_eigenvalues.tsv"))

  pc1 <- round(100 * eig$rel_eig[1], 1)
  pc2 <- round(100 * eig$rel_eig[2], 1)

  plot_scatter_ellipse(
    coords, "Axis.1", "Axis.2",
    paste0(title_prefix, " PCoA (", corr, "; neg_eig_raw=", neg, ")"),
    paste0("PCoA1 (", pc1, "%)"),
    paste0("PCoA2 (", pc2, "%)"),
    file.path(out_dir, "figures", "PCoA_PC1_PC2.pdf"),
    file.path(out_dir, "figures", "PCoA_PC1_PC2.png")
  )

  set.seed(1)
  nmds <- metaMDS(dist_obj, k = 2, trymax = 200, trace = 0)

  nsc_raw <- as.data.frame(scores(nmds, display = "sites")) |>
    rownames_to_column("sample_id")
  nsc <- attach_meta(nsc_raw, meta_map, out_dir = out_dir, tag = "nmds")

  write_tsv0(nsc, file.path(out_dir, "tables", "nmds_scores_per_sample.tsv"))
  write_tsv0(tibble(stress = nmds$stress, converged = nmds$converged, tries = nmds$tries),
             file.path(out_dir, "tables", "nmds_stress.tsv"))

  plot_scatter_ellipse(
    nsc, "NMDS1", "NMDS2",
    paste0(title_prefix, " NMDS (stress=", round(nmds$stress, 4), ")"),
    "NMDS1", "NMDS2",
    file.path(out_dir, "figures", "NMDS_1_2.pdf"),
    file.path(out_dir, "figures", "NMDS_1_2.png")
  )
}

clr_transform <- function(mat_feat_by_samp, pseudocount = 0.5) {
  x <- mat_feat_by_samp + pseudocount
  x <- sweep(x, 2, colSums(x), "/")
  lx <- log(x)
  sweep(lx, 2, colMeans(lx), "-")
}

## ---------- output root ----------
out_root <- "results_16S_ASV_beta"
dir.create(out_root, showWarnings = FALSE)
dir.create(file.path(out_root, "00_manifest"), showWarnings = FALSE)

## ---------- locate inputs ----------
meta_file <- detect_file("metadata.csv", c("^metadata\\.csv$"))
asv_file  <- detect_file("table_f_even_depth_w_tax.xls", c("table_f_even_depth_w_tax\\.(xls|tsv|txt)$"))
tree_file <- detect_file(NULL, c("rooted.*\\.(tre|nwk)$", ".*tree.*\\.(tre|nwk)$", ".*\\.(tre|nwk)$"))

## ---------- load metadata ----------
meta_raw <- read_csv_smart(meta_file) |>
  janitor::clean_names(replace = c("\u03bc" = "u"))

sample_col <- find_col(meta_raw, c("^sample$"))
group_col  <- find_col(meta_raw, c("^group$"))
age_col    <- find_col(meta_raw, c("^age$"))
bmi_col    <- find_col(meta_raw, c("^bmi"))
dur_col    <- find_col(meta_raw, c("disease_duration"))
con_col    <- find_col(meta_raw, c("constipation"))
ldopa_col  <- find_col(meta_raw, c("levodopa_daily_dose"))

meta <- meta_raw |>
  mutate(sample_id = normalize_id(.data[[sample_col]]),
         group = .data[[group_col]]) |>
  filter(!is.na(sample_id) & sample_id != "") |>
  filter(group %in% GROUP_LEVELS) |>
  transmute(
    sample_id,
    group = factor(group, levels = GROUP_LEVELS),
    age = as.numeric(.data[[age_col]]),
    bmi = as.numeric(.data[[bmi_col]]),
    disease_duration = as.numeric(.data[[dur_col]]),
    constipation = factor(as.integer(.data[[con_col]]), levels = c(0,1), labels = c("No","Yes")),
    levodopa_dose = if (!is.na(ldopa_col)) as.numeric(.data[[ldopa_col]]) else NA_real_
  )
stopifnot(nrow(meta) == 40)
meta_rn <- meta |> column_to_rownames("sample_id")
write_tsv0(meta, file.path(out_root, "00_manifest", "metadata_model_40.tsv"))

## ---------- read 16S table ----------
hdr <- fread(asv_file, sep = "\t", nrows = 0, header = TRUE, data.table = FALSE, check.names = FALSE)
cols <- names(hdr)
id_col <- if ("#ID" %in% cols) "#ID" else cols[1]
tax_col <- if ("taxonomy" %in% cols) "taxonomy" else NA_character_

sample_cols <- cols[str_detect(cols, "^PD\\d+$")]
keep_samples <- intersect(meta$sample_id, sample_cols)
stopifnot(length(keep_samples) == 40)
keep_samples <- meta$sample_id

sel <- c(id_col, keep_samples, if (!is.na(tax_col)) tax_col)
asv_dt <- fread(asv_file, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE, select = sel)

asv_id <- as.character(asv_dt[[id_col]])
mat_counts <- asv_dt[, keep_samples, drop = FALSE]
mat_counts[] <- lapply(mat_counts, as.numeric)
mat_counts <- as.matrix(mat_counts)
rownames(mat_counts) <- asv_id

if (!is.na(tax_col)) {
  write_tsv0(tibble(asv_id = asv_id, taxonomy = asv_dt[[tax_col]]),
             file.path(out_root, "00_manifest", "ASV_taxonomy_table.tsv"))
}

write_tsv0(tibble(file = c(meta_file, asv_file, tree_file),
                  role = c("metadata", "16S_ASV_counts", "rooted_tree")),
           file.path(out_root, "00_manifest", "inputs_manifest.tsv"))

mat_counts_nz <- mat_counts[rowSums(mat_counts) > 0, , drop = FALSE]
write_tsv0(tibble(asv_total = nrow(mat_counts),
                  asv_nonzero = nrow(mat_counts_nz),
                  samples = ncol(mat_counts_nz)),
           file.path(out_root, "00_manifest", "ASV_nonzero_summary.tsv"))

## ---------- matrices for distances ----------
mat_rel <- sweep(mat_counts_nz, 2, colSums(mat_counts_nz), "/")
mat_pa  <- (mat_counts_nz > 0) * 1
mat_clr <- clr_transform(mat_counts_nz, pseudocount = 0.5)

saveRDS(mat_rel, file.path(out_root, "00_manifest", "mat_rel_TSS.rds"), compress = "xz")
saveRDS(mat_pa,  file.path(out_root, "00_manifest", "mat_pa.rds"), compress = "xz")
saveRDS(mat_clr, file.path(out_root, "00_manifest", "mat_clr.rds"), compress = "xz")

run_distance_branch <- function(dist_name, dist_obj, title_prefix) {
  out_dir <- file.path(out_root, dist_name)
  dir_init(out_dir)

  write_tsv0(tibble(distance = dist_name, samples = nrow(meta_rn), permutations = 999),
             file.path(out_dir, "00_manifest", "params.tsv"))
  write_tsv0(meta |> select(sample_id, group), file.path(out_dir, "00_manifest", "sample_list.tsv"))

  write_dist_tables(dist_obj, file.path(out_dir, "tables", "dist"))
  run_ordination_outputs(dist_obj, meta_rn, out_dir, title_prefix)
  run_adonis_suite(dist_obj, meta_rn, out_dir, perm = 999)
  run_betadisper_suite(dist_obj, meta_rn, out_dir, perm = 999)

  capture.output(sessionInfo(), file = file.path(out_dir, "logs", "sessionInfo.txt"))
}

run_distance_branch("Bray_Curtis", vegdist(t(mat_rel), method = "bray"), "16S ASV")
run_distance_branch("Jaccard_presence_absence", vegdist(t(mat_pa), method = "jaccard", binary = TRUE), "16S ASV")
run_distance_branch("Aitchison_CLR_Euclidean", dist(t(mat_clr), method = "euclidean"), "16S ASV")

## UniFrac
tree <- read.tree(tree_file)
tree$tip.label <- gsub("^'+|'+$", "", tree$tip.label)

keep_taxa <- rownames(mat_counts_nz)
tree <- drop.tip(tree, setdiff(tree$tip.label, keep_taxa))
tree <- multi2di(tree, random = FALSE)

ps <- phyloseq(otu_table(mat_counts_nz, taxa_are_rows = TRUE),
               sample_data(meta_rn),
               phy_tree(tree))
ps <- prune_taxa(taxa_names(ps) %in% tree$tip.label, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

run_distance_branch("UniFrac_unweighted", UniFrac(ps, weighted = FALSE, normalized = TRUE, parallel = FALSE), "16S ASV")
run_distance_branch("UniFrac_weighted",   UniFrac(ps, weighted = TRUE,  normalized = TRUE, parallel = FALSE), "16S ASV")

message("16S beta full suite done -> ", out_root)