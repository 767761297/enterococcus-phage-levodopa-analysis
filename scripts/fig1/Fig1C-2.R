## =========================
## 06_heatmap_16S_Genus_p05_3models.R
## Genus-level heatmap (3 models): CLR (or residual CLR) -> row Z-score -> heatmap
## Inclusion rule: p < 0.05 (Wilcoxon), BH computed but NOT used for inclusion here.
## Inputs in working dir:
##   - metadata.csv
##   - table_filtered_w_tax.xls 或  - table_f_even_depth_w_tax.xls
## Output:
##   - results_16S_Genus_heatmap_p05/ 或 - results_16S_Genus_heatmap_p05_even/
## =========================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(readr)
  library(janitor)
  library(tibble)
  library(ggplot2)
})

## ---------- user params ----------
META_FILE <- "metadata.csv"
INPUT_16S <- "table_f_even_depth_w_tax.xls" #抽平/未抽平

OUT_ROOT <- "results_16S_Genus_heatmap_p05_even" #抽平/未抽平

GROUP_LEVELS <- c("PD_GR", "PD_MR")  # GR first, MR second

# genus aggregation/filter (match your DA defaults; you can tighten if needed)
MIN_PREVALENCE_ALL <- 0.10
MIN_TOTAL_COUNT <- 10

# CLR
PSEUDOCOUNT <- 0.5

# heatmap inclusion
P_CUTOFF <- 0.05
FALLBACK_TOP_N <- 30  # if no p < 0.05 taxa

# plot
CLUSTER_ROWS <- TRUE
CLUSTER_COLS <- FALSE   # keep sample order by group then sample_id
MAKE_PDF <- TRUE
MAKE_PNG <- TRUE

## ---------- helpers ----------
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
fill_unknown <- function(x, label) {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "NA"] <- paste0("Unassigned_", label)
  x
}
parse_tax_16s <- function(tax) {
  parts <- str_split(as.character(tax), ";", simplify = TRUE)
  get_rank <- function(row, key) {
    cand <- row[str_detect(row, paste0("^\\s*", key, "_+"))]
    if (length(cand) == 0) return(NA_character_)
    v <- str_trim(cand[1])
    v <- str_replace(v, paste0("^\\s*", key, "_+"), "")
    if (v == "" || v == "NA") NA_character_ else v
  }
  tibble(
    domain  = apply(parts, 1, get_rank, key = "d"),
    phylum  = apply(parts, 1, get_rank, key = "p"),
    class   = apply(parts, 1, get_rank, key = "c"),
    order   = apply(parts, 1, get_rank, key = "o"),
    family  = apply(parts, 1, get_rank, key = "f"),
    genus   = apply(parts, 1, get_rank, key = "g"),
    species = apply(parts, 1, get_rank, key = "s")
  )
}
aggregate_counts <- function(mat_feat_by_sample, feature_to_group_named) {
  grp <- feature_to_group_named[rownames(mat_feat_by_sample)]
  rowsum(mat_feat_by_sample, group = grp, reorder = TRUE)
}
clr_transform <- function(mat_feat_by_sample, pseudocount = 0.5) {
  x <- mat_feat_by_sample + pseudocount
  x <- sweep(x, 2, colSums(x), "/")
  lx <- log(x)
  sweep(lx, 2, colMeans(lx), "-")
}
residualize_matrix <- function(Y_samp_by_feat, X_design) {
  qrX <- qr(X_design)
  qr.resid(qrX, Y_samp_by_feat)
}
row_zscore <- function(mat_feat_by_sample) {
  # features x samples
  mu <- rowMeans(mat_feat_by_sample, na.rm = TRUE)
  sdv <- apply(mat_feat_by_sample, 1, sd, na.rm = TRUE)
  sdv[sdv == 0 | is.na(sdv)] <- 1
  sweep(sweep(mat_feat_by_sample, 1, mu, "-"), 1, sdv, "/")
}

# robust Wilcoxon: compute U from ranks (rank-biserial in [-1,1])
wilcox_one <- function(x_mr, x_gr) {
  n1 <- length(x_mr); n2 <- length(x_gr)
  wt <- suppressWarnings(wilcox.test(x_mr, x_gr, exact = FALSE))
  p <- wt$p.value
  Wstat <- as.numeric(wt$statistic)

  z <- c(x_mr, x_gr)
  r <- rank(z, ties.method = "average")
  Wrank <- sum(r[seq_len(n1)])
  U <- Wrank - n1 * (n1 + 1) / 2
  rb <- (2 * U) / (n1 * n2) - 1
  rb <- max(-1, min(1, rb))
  med_diff <- median(x_mr, na.rm = TRUE) - median(x_gr, na.rm = TRUE)

  c(p_value = p, W = Wstat, U = U, rank_biserial = rb, median_diff = med_diff)
}

plot_heatmap_pheatmap <- function(mat_z, anno_col, out_pdf, out_png, cluster_rows, cluster_cols) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) return(FALSE)

  # pheatmap wants annotation as data.frame with rownames=colnames(mat)
  ann <- as.data.frame(anno_col)
  rownames(ann) <- colnames(mat_z)

  if (MAKE_PDF) {
    pdf(out_pdf, width = 9, height = max(4, 0.18 * nrow(mat_z) + 2))
    pheatmap::pheatmap(mat_z,
                       cluster_rows = cluster_rows,
                       cluster_cols = cluster_cols,
                       annotation_col = ann,
                       show_colnames = FALSE,
                       fontsize_row = 8)
    dev.off()
  }
  if (MAKE_PNG) {
    png(out_png, width = 2200, height = max(1400, 40 * nrow(mat_z) + 600), res = 300)
    pheatmap::pheatmap(mat_z,
                       cluster_rows = cluster_rows,
                       cluster_cols = cluster_cols,
                       annotation_col = ann,
                       show_colnames = FALSE,
                       fontsize_row = 8)
    dev.off()
  }
  TRUE
}

plot_heatmap_gg <- function(mat_z, sample_order, out_pdf, out_png) {
  df <- as.data.frame(mat_z) |>
    rownames_to_column("feature") |>
    pivot_longer(-feature, names_to = "sample_id", values_to = "z") |>
    mutate(sample_id = factor(sample_id, levels = sample_order))

  # order features by clustering (optional)
  feat_order <- rownames(mat_z)
  if (CLUSTER_ROWS && nrow(mat_z) >= 2) {
    hc <- hclust(dist(mat_z))
    feat_order <- rownames(mat_z)[hc$order]
  }
  df$feature <- factor(df$feature, levels = rev(feat_order))

  p <- ggplot(df, aes(x = sample_id, y = feature, fill = z)) +
    geom_tile() +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    labs(x = NULL, y = NULL, fill = "Z")

  if (MAKE_PDF) ggsave(out_pdf, p, width = 9, height = max(4, 0.18 * nlevels(df$feature) + 2))
  if (MAKE_PNG) ggsave(out_png, p, width = 9, height = max(4, 0.18 * nlevels(df$feature) + 2), dpi = 300)
}

## ---------- load metadata ----------
meta_raw <- read_csv_smart(META_FILE) |> janitor::clean_names(replace = c("\u03bc" = "u"))

sample_col <- find_col(meta_raw, c("^sample$"))
group_col  <- find_col(meta_raw, c("^group$"))
age_col    <- find_col(meta_raw, c("^age$"))
bmi_col    <- find_col(meta_raw, c("^bmi"))
dur_col    <- find_col(meta_raw, c("disease_duration"))
con_col    <- find_col(meta_raw, c("constipation"))
dose_col   <- find_col(meta_raw, c("levodopa.*dose"))

stopifnot(!is.na(sample_col), !is.na(group_col))

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
    levodopa_dose = if (!is.na(dose_col)) as.numeric(.data[[dose_col]]) else NA_real_
  )

stopifnot(nrow(meta) == 40)

# sample order: GR then MR, within group by sample_id
meta <- meta |>
  arrange(group, sample_id)
sample_order <- meta$sample_id
meta_rn <- meta |> column_to_rownames("sample_id")

## ---------- read 16S table (tab-delimited with .xls extension) ----------
dt <- fread(INPUT_16S, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)

id_col  <- if ("#ID" %in% names(dt)) "#ID" else names(dt)[1]
tax_col <- if ("taxonomy" %in% names(dt)) "taxonomy" else if ("taxon" %in% names(dt)) "taxon" else NA_character_
stopifnot(!is.na(tax_col))

sample_cols <- names(dt)[str_detect(names(dt), "^PD\\d+$")]
keep_samples <- intersect(sample_order, sample_cols)
stopifnot(length(keep_samples) == 40)
keep_samples <- sample_order

asv_id <- as.character(dt[[id_col]])
mat_asv <- dt[, keep_samples, drop = FALSE]
mat_asv[] <- lapply(mat_asv, as.numeric)
mat_asv <- as.matrix(mat_asv)
rownames(mat_asv) <- asv_id

tax_tbl <- tibble(asv_id = asv_id, taxonomy_raw = dt[[tax_col]]) |>
  bind_cols(parse_tax_16s(dt[[tax_col]])) |>
  mutate(genus = fill_unknown(genus, "genus"))

## ---------- aggregate to genus & filter ----------
map_asv_genus <- setNames(tax_tbl$genus, tax_tbl$asv_id)
mat_genus <- aggregate_counts(mat_asv, map_asv_genus)

prev_all <- rowMeans(mat_genus > 0)
tot <- rowSums(mat_genus)
keep <- (prev_all >= MIN_PREVALENCE_ALL) & (tot >= MIN_TOTAL_COUNT)
mat_genus_f <- mat_genus[keep, , drop = FALSE]

dir.create(OUT_ROOT, showWarnings = FALSE)
write_tsv0(tibble(
  n_total = nrow(mat_genus),
  n_keep = nrow(mat_genus_f),
  min_prevalence_all = MIN_PREVALENCE_ALL,
  min_total_count = MIN_TOTAL_COUNT
), file.path(OUT_ROOT, "00_manifest_filter_summary_Genus.tsv"))

saveRDS(mat_genus_f, file.path(OUT_ROOT, "genus_counts_filtered.rds"), compress = "xz")

## ---------- compute CLR ----------
clr_genus <- clr_transform(mat_genus_f, pseudocount = PSEUDOCOUNT)  # features x samples
saveRDS(clr_genus, file.path(OUT_ROOT, "genus_CLR_features_by_sample.rds"), compress = "xz")

## ---------- model definitions ----------
models <- list(
  unadjusted = list(type = "clr", X = NULL),
  primary = list(type = "resid",
                 X = model.matrix(~ age + bmi + disease_duration + constipation, data = meta)),
  sensitivity_add_dose = list(type = "resid",
                              X = model.matrix(~ age + bmi + disease_duration + constipation + levodopa_dose, data = meta))
)

## ---------- run per model ----------
for (mname in names(models)) {
  out_dir <- file.path(OUT_ROOT, mname)
  dir_init(out_dir)
  write_tsv0(tibble(model = mname, p_cutoff = P_CUTOFF,
                    cluster_rows = CLUSTER_ROWS, cluster_cols = CLUSTER_COLS),
             file.path(out_dir, "00_manifest", "params.tsv"))

  # build test matrix (samples x features)
  Y <- t(clr_genus)  # samples x features
  rownames(Y) <- colnames(clr_genus)
  colnames(Y) <- rownames(clr_genus)

  if (models[[mname]]$type == "resid") {
    Y <- residualize_matrix(Y, models[[mname]]$X)
  }

  # compute wilcoxon per feature (MR vs GR)
  idx_gr <- which(meta$group == "PD_GR")
  idx_mr <- which(meta$group == "PD_MR")

  feats <- colnames(Y)
  stats <- t(vapply(feats, function(f) {
    wilcox_one(Y[idx_mr, f], Y[idx_gr, f])
  }, numeric(5)))

  res <- tibble(
    feature = feats,
    effect_median_diff = stats[, "median_diff"],
    W = stats[, "W"],
    U = stats[, "U"],
    rank_biserial = stats[, "rank_biserial"],
    p_value = stats[, "p_value"]
  ) |>
    mutate(q_value = p.adjust(p_value, method = "BH"),
           direction = case_when(effect_median_diff > 0 ~ "PD_MR>PD_GR",
                                 effect_median_diff < 0 ~ "PD_MR<PD_GR",
                                 TRUE ~ "NA"),
           level = "Genus")

  write_tsv0(res, file.path(out_dir, "tables", "wilcoxon_all_results.tsv"))

  # select p<0.05; fallback to top N if none
  sel <- res |> filter(p_value < P_CUTOFF) |> arrange(p_value)
  note <- ""
  if (nrow(sel) == 0) {
    sel <- res |> arrange(p_value) |> slice_head(n = min(FALLBACK_TOP_N, nrow(res)))
    note <- paste0("No features with p<", P_CUTOFF, "; fallback to top ", nrow(sel), " by p.")
  } else {
    note <- paste0("Selected features: p<", P_CUTOFF, " (n=", nrow(sel), ").")
  }
  writeLines(note, con = file.path(out_dir, "logs", "selection_note.txt"))
  write_tsv0(sel, file.path(out_dir, "tables", "selected_by_p.tsv"))

  sel_feats <- sel$feature

  # matrices for heatmap
  mat_test_feat_by_samp <- t(Y[, sel_feats, drop = FALSE])  # features x samples
  mat_test_feat_by_samp <- mat_test_feat_by_samp[, sample_order, drop = FALSE]

  mat_z <- row_zscore(mat_test_feat_by_samp)

  # output reproducible tables
  write_tsv0(as.data.frame(mat_test_feat_by_samp) |> rownames_to_column("feature"),
             file.path(out_dir, "tables", "heatmap_matrix_testvalue_wide.tsv"))
  write_tsv0(as.data.frame(mat_z) |> rownames_to_column("feature"),
             file.path(out_dir, "tables", "heatmap_matrix_zscore_wide.tsv"))
  write_tsv0(meta, file.path(out_dir, "tables", "sample_annotation.tsv"))

  # per-sample long (for reproducing any plot)
  long <- as.data.frame(mat_test_feat_by_samp) |>
    rownames_to_column("feature") |>
    pivot_longer(cols = all_of(sample_order), names_to = "sample_id", values_to = "test_value") |>
    left_join(as.data.frame(mat_z) |> rownames_to_column("feature") |>
                pivot_longer(cols = all_of(sample_order), names_to = "sample_id", values_to = "z_value"),
              by = c("feature","sample_id")) |>
    left_join(meta, by = "sample_id") |>
    left_join(sel |> select(feature, p_value, q_value, effect_median_diff, rank_biserial, direction),
              by = "feature")
  write_tsv0(long, file.path(out_dir, "tables", "per_sample_heatmap_long.tsv"))

  # plot heatmap
  anno_col <- data.frame(group = meta$group)
  rownames(anno_col) <- meta$sample_id
  anno_col <- anno_col[sample_order, , drop = FALSE]

  out_pdf <- file.path(out_dir, "figures", paste0("Genus_heatmap_p", P_CUTOFF, "_", mname, ".pdf"))
  out_png <- file.path(out_dir, "figures", paste0("Genus_heatmap_p", P_CUTOFF, "_", mname, ".png"))

  ok <- plot_heatmap_pheatmap(mat_z, anno_col, out_pdf, out_png, CLUSTER_ROWS, CLUSTER_COLS)
  if (!ok) {
    # fallback to ggplot heatmap
    plot_heatmap_gg(mat_z, sample_order, out_pdf, out_png)
  }

  message("Done heatmap: ", mname, " | features=", nrow(sel))
}

message("All genus heatmaps done -> ", OUT_ROOT)