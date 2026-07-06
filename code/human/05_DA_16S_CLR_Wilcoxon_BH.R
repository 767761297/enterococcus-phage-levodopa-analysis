## =========================
## 05_DA_16S_CLR_Wilcoxon_BH_UPGRADED.R
## 16S DA: CLR + Wilcoxon + BH
## Levels: ASV + Genus + Species
## Models: unadjusted / primary / sensitivity(+Levodopa dose) via residualization
## Upgrades:
##   - robust rank-biserial in [-1,1] via U computed from ranks
##   - optional independent filtering: mean_RA and/or min prevalence within each group
##   - ASV-level presence/absence supplementary test (Fisher/Chi-square + BH)
## Inputs in working dir:
##   - metadata.csv
##   - table_filtered_w_tax.xls   或  - table_f_even_depth_w_tax.xls
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

## ---------- REQUIRED INPUTS ----------
META_FILE <- "metadata.csv"
INPUT_16S <- "table_f_even_depth_w_tax.xls"

## ---------- OUTPUT ----------
OUT_ROOT <- "results_16S_DA_CLR_Wilcoxon_BH_even" 
# 输出文件夹为results_16S_DA_CLR_Wilcoxon_BH_even 或results_16S_DA_CLR_Wilcoxon_BH
#输出文件夹名称根据输入文件名称

## ---------- GLOBAL SETTINGS ----------
GROUP_LEVELS <- c("PD_GR", "PD_MR")   # effect defined as MR - GR
FDR_CUTOFF <- 0.10

# Base filtering (before CLR): applied to each level matrix
MIN_PREVALENCE_ALL <- 0.10     # overall prevalence >= 10% samples nonzero
MIN_TOTAL_COUNT <- 10          # total counts across all samples

# CLR pseudocount for counts
PSEUDOCOUNT <- 0.5

# Export per-sample long data for: sig + top N by q
TOP_N_BY_Q <- 50

# Optional plots
MAKE_VOLCANO_PLOTS <- TRUE

# Optional nominal candidates (exploratory; do NOT treat as main significant)
EXPORT_NOMINAL_CANDIDATES <- TRUE
NOMINAL_P_CUTOFF <- 0.05
MIN_ABS_RB_CANDIDATE <- 0.30

## ---------- OPTIONAL INDEPENDENT FILTERING (OFF by default) ----------
## Set to NA to disable; otherwise enforce extra filtering criteria.
## 1) Require BOTH groups to have prevalence >= threshold
MIN_PREV_BY_GROUP_ASV <- NA_real_
MIN_PREV_BY_GROUP_GENUS <- NA_real_
MIN_PREV_BY_GROUP_SPECIES <- NA_real_

## 2) Require mean relative abundance >= threshold (across all samples)
## (mean of count/total_count_per_sample)
MIN_MEAN_RA_ASV <- NA_real_
MIN_MEAN_RA_GENUS <- NA_real_
MIN_MEAN_RA_SPECIES <- NA_real_

## ---------- ASV PRESENCE/ABSENCE SUPPLEMENT ----------
RUN_ASV_PRESENCE_SUPPLEMENT <- TRUE

# Presence filter (to keep runtime reasonable). Applied on ASV counts (unaggregated).
PRES_MIN_PREVALENCE_ALL <- 0.05      # overall presence >= 5% samples
PRES_MIN_TOTAL_COUNT <- 5            # total counts >= 5

# Method: "auto" (fisher if <=5000 features else chisq), or "fisher", or "chisq"
PRES_TEST_METHOD <- "auto"

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

safe_mode <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != "" & x != "NA"]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

## taxonomy parser (supports d__ / d_ style)
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
  stopifnot(all(rownames(mat_feat_by_sample) %in% names(feature_to_group_named)))
  grp <- feature_to_group_named[rownames(mat_feat_by_sample)]
  rowsum(mat_feat_by_sample, group = grp, reorder = TRUE)
}

# Efficient mean RA: mean_j (x_ij / total_j) = (mat %*% (1/total)) / n
mean_ra_fast <- function(mat_feat_by_sample) {
  tot <- colSums(mat_feat_by_sample)
  w <- 1 / tot
  as.numeric(mat_feat_by_sample %*% w) / length(w)
}

filter_features <- function(mat_feat_by_sample, meta,
                            min_prev_all, min_total_count,
                            min_prev_by_group = NA_real_,
                            min_mean_ra = NA_real_) {
  # overall prevalence and total
  prev_all <- rowMeans(mat_feat_by_sample > 0)
  tot <- rowSums(mat_feat_by_sample)

  keep <- (prev_all >= min_prev_all) & (tot >= min_total_count)

  # group prevalence (optional)
  if (!is.na(min_prev_by_group)) {
    idx_gr <- which(meta$group == "PD_GR")
    idx_mr <- which(meta$group == "PD_MR")
    prev_gr <- rowMeans(mat_feat_by_sample[, idx_gr, drop = FALSE] > 0)
    prev_mr <- rowMeans(mat_feat_by_sample[, idx_mr, drop = FALSE] > 0)
    keep <- keep & (prev_gr >= min_prev_by_group) & (prev_mr >= min_prev_by_group)
  } else {
    prev_gr <- prev_mr <- rep(NA_real_, nrow(mat_feat_by_sample))
  }

  # mean RA (optional)
  if (!is.na(min_mean_ra)) {
    mra <- mean_ra_fast(mat_feat_by_sample)
    keep <- keep & (mra >= min_mean_ra)
  } else {
    mra <- rep(NA_real_, nrow(mat_feat_by_sample))
  }

  summary <- tibble(
    n_total = nrow(mat_feat_by_sample),
    n_keep  = sum(keep),
    min_prev_all = min_prev_all,
    min_total_count = min_total_count,
    min_prev_by_group = min_prev_by_group,
    min_mean_ra = min_mean_ra
  )

  list(
    mat = mat_feat_by_sample[keep, , drop = FALSE],
    summary = summary,
    diagnostics = tibble(
      feature = rownames(mat_feat_by_sample),
      prev_all = prev_all,
      prev_gr = prev_gr,
      prev_mr = prev_mr,
      total_count = tot,
      mean_ra = mra,
      keep = keep
    )
  )
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

# Robust Wilcoxon summary: compute U from ranks (ties handled) => rank-biserial in [-1,1]
wilcox_one <- function(x_mr, x_gr) {
  n1 <- length(x_mr); n2 <- length(x_gr)

  wt <- suppressWarnings(wilcox.test(x_mr, x_gr, exact = FALSE))
  W_stat <- as.numeric(wt$statistic)
  p <- wt$p.value

  # U from rank sum (independent of what wt$statistic encodes)
  z <- c(x_mr, x_gr)
  r <- rank(z, ties.method = "average")
  W_rank_sum <- sum(r[seq_len(n1)])
  U <- W_rank_sum - n1 * (n1 + 1) / 2

  rb <- (2 * U) / (n1 * n2) - 1
  rb <- max(-1, min(1, rb))  # clamp numeric jitter
  med_diff <- median(x_mr, na.rm = TRUE) - median(x_gr, na.rm = TRUE)

  c(p_value = p, W = W_stat, U = U, rank_biserial = rb, median_diff = med_diff)
}

run_wilcox_suite <- function(level_name, mat_counts, anno_tbl, meta, model_name) {
  out_dir <- file.path(OUT_ROOT, level_name, model_name)
  dir_init(out_dir)

  # CLR
  clr_mat <- clr_transform(mat_counts, pseudocount = PSEUDOCOUNT)   # features x samples
  saveRDS(clr_mat, file.path(out_dir, "00_manifest", "clr_matrix_features_by_sample.rds"), compress = "xz")

  # testing matrix: samples x features
  Y <- t(clr_mat)
  rownames(Y) <- colnames(mat_counts)
  colnames(Y) <- rownames(mat_counts)

  # residualization for adjusted models
  if (model_name == "primary") {
    X <- model.matrix(~ age + bmi + disease_duration + constipation, data = meta)
    Y <- residualize_matrix(Y, X)
  }
  if (model_name == "sensitivity_add_dose") {
    X <- model.matrix(~ age + bmi + disease_duration + constipation + levodopa_dose, data = meta)
    Y <- residualize_matrix(Y, X)
  }

  saveRDS(Y, file.path(out_dir, "00_manifest", "test_matrix_samples_by_feature.rds"), compress = "xz")

  grp <- meta$group
  idx_gr <- which(grp == "PD_GR")
  idx_mr <- which(grp == "PD_MR")

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
    mutate(
      q_value = p.adjust(p_value, method = "BH"),
      direction = case_when(effect_median_diff > 0 ~ "PD_MR>PD_GR",
                            effect_median_diff < 0 ~ "PD_MR<PD_GR",
                            TRUE ~ "NA")
    ) |>
    left_join(anno_tbl, by = "feature")

  write_tsv0(res, file.path(out_dir, "tables", "all_results.tsv"))

  sig <- res |> filter(!is.na(q_value) & q_value < FDR_CUTOFF)
  write_tsv0(sig, file.path(out_dir, "tables", paste0("significant_FDR_", FDR_CUTOFF, ".tsv")))

  # exploratory candidates
  if (EXPORT_NOMINAL_CANDIDATES) {
    cand <- res |>
      filter(!is.na(p_value) & p_value < NOMINAL_P_CUTOFF) |>
      filter(abs(rank_biserial) >= MIN_ABS_RB_CANDIDATE) |>
      arrange(p_value)
    write_tsv0(cand, file.path(out_dir, "tables", paste0("candidates_p_lt_", NOMINAL_P_CUTOFF,
                                                        "_absRB_ge_", MIN_ABS_RB_CANDIDATE, ".tsv")))
  }

  # volcano data
  volc <- res |>
    transmute(feature,
              effect = effect_median_diff,
              neglog10_q = -log10(pmax(q_value, 1e-300)),
              q_value,
              direction,
              label = ifelse(q_value < FDR_CUTOFF, feature, ""))
  write_tsv0(volc, file.path(out_dir, "tables", "volcano_data.tsv"))

  # select features for per-sample export: sig + top N
  top_feats <- res |> arrange(q_value) |> slice_head(n = TOP_N_BY_Q) |> pull(feature)
  sig_feats <- sig |> pull(feature)
  selected <- unique(c(sig_feats, top_feats))
  selected <- selected[selected %in% rownames(mat_counts)]
  write_tsv0(tibble(feature = selected), file.path(out_dir, "tables", "selected_features.tsv"))

  # per-sample long (counts + RA + CLR + test_value)
  if (length(selected) > 0) {
    cnt_long <- as.data.frame(mat_counts[selected, , drop = FALSE]) |>
      rownames_to_column("feature") |>
      pivot_longer(cols = all_of(colnames(mat_counts)),
                   names_to = "sample_id", values_to = "count")

    samp_tot <- colSums(mat_counts)
    ra_long <- cnt_long |>
      mutate(relative_abundance = count / samp_tot[sample_id])

    clr_sel <- as.data.frame(clr_mat[selected, , drop = FALSE]) |>
      rownames_to_column("feature") |>
      pivot_longer(cols = all_of(colnames(mat_counts)),
                   names_to = "sample_id", values_to = "clr_value")

    test_sel <- as.data.frame(t(Y[, selected, drop = FALSE])) |>
      rownames_to_column("feature")
    colnames(test_sel)[-1] <- rownames(Y)
    test_long <- test_sel |>
      pivot_longer(cols = all_of(rownames(Y)),
                   names_to = "sample_id", values_to = "test_value")

    out_long <- ra_long |>
      left_join(clr_sel, by = c("feature","sample_id")) |>
      left_join(test_long, by = c("feature","sample_id")) |>
      left_join(meta |> select(sample_id, group, age, bmi, disease_duration, constipation, levodopa_dose),
                by = "sample_id") |>
      left_join(res |> select(feature, effect_median_diff, rank_biserial, p_value, q_value, direction),
                by = "feature")

    write_tsv0(out_long, file.path(out_dir, "tables", "per_sample_selected_long.tsv"))
    write_tsv0(test_sel, file.path(out_dir, "tables", "heatmap_matrix_selected_wide.tsv"))
  }

  if (MAKE_VOLCANO_PLOTS) {
    p <- ggplot(volc, aes(x = effect, y = neglog10_q)) +
      geom_point(alpha = 0.7) +
      theme_bw() +
      geom_vline(xintercept = 0, linetype = 2) +
      labs(title = paste0(level_name, " ", model_name, " volcano"),
           x = "Median difference (MR - GR) on test scale",
           y = "-log10(q)")
    ggsave(file.path(out_dir, "figures", "volcano.png"), p, width = 6.5, height = 5, dpi = 300)
    ggsave(file.path(out_dir, "figures", "volcano.pdf"), p, width = 6.5, height = 5)
  }

  capture.output(sessionInfo(), file = file.path(out_dir, "logs", "sessionInfo.txt"))
}

## ---------- ASV presence/absence supplement ----------
run_asv_presence_supplement <- function(mat_asv, meta, tax_tbl) {
  out_dir <- file.path(OUT_ROOT, "ASV", "presence_supplement")
  dir_init(out_dir)

  # filter for feasibility
  prev_all <- rowMeans(mat_asv > 0)
  tot <- rowSums(mat_asv)
  keep <- (prev_all >= PRES_MIN_PREVALENCE_ALL) & (tot >= PRES_MIN_TOTAL_COUNT)
  mat_p <- mat_asv[keep, , drop = FALSE]

  write_tsv0(tibble(
    n_total = nrow(mat_asv),
    n_keep  = nrow(mat_p),
    pres_min_prev_all = PRES_MIN_PREVALENCE_ALL,
    pres_min_total_count = PRES_MIN_TOTAL_COUNT,
    pres_test_method = PRES_TEST_METHOD
  ), file.path(out_dir, "00_manifest", "presence_filter_params.tsv"))

  if (nrow(mat_p) == 0) {
    writeLines("No ASVs passed presence filter; presence supplement skipped.",
               con = file.path(out_dir, "logs", "note.txt"))
    return(invisible(NULL))
  }

  idx_gr <- which(meta$group == "PD_GR")
  idx_mr <- which(meta$group == "PD_MR")
  n_gr <- length(idx_gr); n_mr <- length(idx_mr)

  pres <- (mat_p > 0) * 1
  a <- rowSums(pres[, idx_mr, drop = FALSE]) # MR present
  c <- rowSums(pres[, idx_gr, drop = FALSE]) # GR present
  b <- n_mr - a
  d <- n_gr - c

  # odds ratio with Haldane-Anscombe correction
  OR <- ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
  logOR <- log(OR)

  # choose test
  method_use <- PRES_TEST_METHOD
  if (method_use == "auto") {
    method_use <- if (nrow(mat_p) <= 5000) "fisher" else "chisq"
  }

  if (method_use == "chisq") {
    n <- n_mr + n_gr
    denom <- (a + b) * (c + d) * (a + c) * (b + d)
    chi2 <- (n * (a * d - b * c)^2) / pmax(denom, 1e-12)
    p <- pchisq(chi2, df = 1, lower.tail = FALSE)
  } else {
    # fisher per feature (can be slow if too many; auto avoids huge n)
    p <- vapply(seq_len(nrow(mat_p)), function(i) {
      m <- matrix(c(a[i], b[i], c[i], d[i]), nrow = 2, byrow = TRUE)
      suppressWarnings(fisher.test(m)$p.value)
    }, numeric(1))
  }

  q <- p.adjust(p, method = "BH")

  cont <- tibble(
    feature = rownames(mat_p),
    MR_present = a, MR_absent = b,
    GR_present = c, GR_absent = d,
    prev_MR = a / n_mr,
    prev_GR = c / n_gr,
    logOR = logOR,
    OR = OR,
    p_value = p,
    q_value = q,
    direction = case_when(logOR > 0 ~ "PD_MR>PD_GR",
                          logOR < 0 ~ "PD_MR<PD_GR",
                          TRUE ~ "NA")
  ) |>
    left_join(tax_tbl |> transmute(feature = asv_id, taxonomy_raw, domain, phylum, class, order, family, genus, species),
              by = "feature")

  write_tsv0(cont, file.path(out_dir, "tables", "asv_presence_all_results.tsv"))
  write_tsv0(cont |> filter(q_value < FDR_CUTOFF),
             file.path(out_dir, "tables", paste0("asv_presence_significant_FDR_", FDR_CUTOFF, ".tsv")))

  if (EXPORT_NOMINAL_CANDIDATES) {
    cand <- cont |> filter(p_value < NOMINAL_P_CUTOFF) |> arrange(p_value)
    write_tsv0(cand, file.path(out_dir, "tables", paste0("asv_presence_candidates_p_lt_", NOMINAL_P_CUTOFF, ".tsv")))
  }

  # per-sample long for top/sig presence candidates
  top_ids <- cont |> arrange(q_value) |> slice_head(n = TOP_N_BY_Q) |> pull(feature)
  sig_ids <- cont |> filter(q_value < FDR_CUTOFF) |> pull(feature)
  pick <- unique(c(top_ids, sig_ids))
  pick <- pick[pick %in% rownames(mat_p)]

  if (length(pick) > 0) {
    long <- as.data.frame(mat_p[pick, , drop = FALSE]) |>
      rownames_to_column("feature") |>
      pivot_longer(cols = colnames(mat_p), names_to = "sample_id", values_to = "count") |>
      mutate(present = as.integer(count > 0)) |>
      left_join(meta, by = "sample_id")
    write_tsv0(long, file.path(out_dir, "tables", "per_sample_presence_top_or_sig_long.tsv"))
  }

  capture.output(sessionInfo(), file = file.path(out_dir, "logs", "sessionInfo.txt"))
}

## ---------- main ----------
dir.create(OUT_ROOT, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "00_manifest"), showWarnings = FALSE)

# metadata
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
meta <- meta |> arrange(sample_id)
write_tsv0(meta, file.path(OUT_ROOT, "00_manifest", "metadata_model_40.tsv"))

write_tsv0(tibble(
  FDR_CUTOFF = FDR_CUTOFF,
  MIN_PREVALENCE_ALL = MIN_PREVALENCE_ALL,
  MIN_TOTAL_COUNT = MIN_TOTAL_COUNT,
  PSEUDOCOUNT = PSEUDOCOUNT,
  TOP_N_BY_Q = TOP_N_BY_Q,
  MIN_PREV_BY_GROUP_ASV = MIN_PREV_BY_GROUP_ASV,
  MIN_PREV_BY_GROUP_GENUS = MIN_PREV_BY_GROUP_GENUS,
  MIN_PREV_BY_GROUP_SPECIES = MIN_PREV_BY_GROUP_SPECIES,
  MIN_MEAN_RA_ASV = MIN_MEAN_RA_ASV,
  MIN_MEAN_RA_GENUS = MIN_MEAN_RA_GENUS,
  MIN_MEAN_RA_SPECIES = MIN_MEAN_RA_SPECIES
), file.path(OUT_ROOT, "00_manifest", "params.tsv"))

# read 16S table
dt <- fread(INPUT_16S, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)

id_col  <- if ("#ID" %in% names(dt)) "#ID" else names(dt)[1]
tax_col <- if ("taxonomy" %in% names(dt)) "taxonomy" else if ("taxon" %in% names(dt)) "taxon" else NA_character_
stopifnot(!is.na(tax_col))

sample_cols <- names(dt)[str_detect(names(dt), "^PD\\d+$")]
keep_samples <- intersect(meta$sample_id, sample_cols)
stopifnot(length(keep_samples) == 40)
keep_samples <- meta$sample_id

asv_id <- as.character(dt[[id_col]])
mat_asv <- dt[, keep_samples, drop = FALSE]
mat_asv[] <- lapply(mat_asv, as.numeric)
mat_asv <- as.matrix(mat_asv)
rownames(mat_asv) <- asv_id

# taxonomy
tax_tbl <- tibble(asv_id = asv_id, taxonomy_raw = dt[[tax_col]]) |>
  bind_cols(parse_tax_16s(dt[[tax_col]])) |>
  mutate(
    genus = fill_unknown(genus, "genus"),
    species = fill_unknown(species, "species"),
    domain = fill_unknown(domain, "domain"),
    phylum = fill_unknown(phylum, "phylum"),
    class  = fill_unknown(class, "class"),
    order  = fill_unknown(order, "order"),
    family = fill_unknown(family, "family")
  )
write_tsv0(tax_tbl, file.path(OUT_ROOT, "00_manifest", "ASV_taxonomy_parsed.tsv"))

## ---------- build levels with upgraded filtering ----------
# ASV
f_asv <- filter_features(mat_asv, meta,
                         min_prev_all = MIN_PREVALENCE_ALL,
                         min_total_count = MIN_TOTAL_COUNT,
                         min_prev_by_group = MIN_PREV_BY_GROUP_ASV,
                         min_mean_ra = MIN_MEAN_RA_ASV)
mat_asv_f <- f_asv$mat
write_tsv0(f_asv$summary |> mutate(level = "ASV"),
           file.path(OUT_ROOT, "00_manifest", "filter_summary_ASV.tsv"))
saveRDS(mat_asv_f, file.path(OUT_ROOT, "00_manifest", "counts_ASV_filtered.rds"), compress = "xz")

# Genus aggregate
map_asv_genus <- setNames(tax_tbl$genus, tax_tbl$asv_id)
mat_genus <- aggregate_counts(mat_asv, map_asv_genus)

f_genus <- filter_features(mat_genus, meta,
                           min_prev_all = MIN_PREVALENCE_ALL,
                           min_total_count = MIN_TOTAL_COUNT,
                           min_prev_by_group = MIN_PREV_BY_GROUP_GENUS,
                           min_mean_ra = MIN_MEAN_RA_GENUS)
mat_genus_f <- f_genus$mat
write_tsv0(f_genus$summary |> mutate(level = "Genus"),
           file.path(OUT_ROOT, "00_manifest", "filter_summary_Genus.tsv"))
saveRDS(mat_genus_f, file.path(OUT_ROOT, "00_manifest", "counts_Genus_filtered.rds"), compress = "xz")
write_tsv0(tibble(asv_id = names(map_asv_genus), genus = unname(map_asv_genus)),
           file.path(OUT_ROOT, "00_manifest", "map_ASV_to_Genus.tsv"))

# Species aggregate
map_asv_species <- setNames(tax_tbl$species, tax_tbl$asv_id)
mat_species <- aggregate_counts(mat_asv, map_asv_species)

f_species <- filter_features(mat_species, meta,
                             min_prev_all = MIN_PREVALENCE_ALL,
                             min_total_count = MIN_TOTAL_COUNT,
                             min_prev_by_group = MIN_PREV_BY_GROUP_SPECIES,
                             min_mean_ra = MIN_MEAN_RA_SPECIES)
mat_species_f <- f_species$mat
write_tsv0(f_species$summary |> mutate(level = "Species"),
           file.path(OUT_ROOT, "00_manifest", "filter_summary_Species.tsv"))
saveRDS(mat_species_f, file.path(OUT_ROOT, "00_manifest", "counts_Species_filtered.rds"), compress = "xz")
write_tsv0(tibble(asv_id = names(map_asv_species), species = unname(map_asv_species)),
           file.path(OUT_ROOT, "00_manifest", "map_ASV_to_Species.tsv"))

## ---------- annotation tables ----------
anno_asv <- tax_tbl |>
  filter(asv_id %in% rownames(mat_asv_f)) |>
  transmute(feature = asv_id, taxonomy_raw, domain, phylum, class, order, family, genus, species)

anno_genus <- tax_tbl |>
  mutate(genus_grp = map_asv_genus[asv_id]) |>
  group_by(genus_grp) |>
  summarise(
    domain = safe_mode(domain),
    phylum = safe_mode(phylum),
    class  = safe_mode(class),
    order  = safe_mode(order),
    family = safe_mode(family),
    .groups = "drop"
  ) |>
  mutate(
    feature = genus_grp,
    genus = genus_grp,
    domain = fill_unknown(domain, "domain"),
    phylum = fill_unknown(phylum, "phylum"),
    class  = fill_unknown(class, "class"),
    order  = fill_unknown(order, "order"),
    family = fill_unknown(family, "family")
  ) |>
  select(feature, domain, phylum, class, order, family, genus) |>
  filter(feature %in% rownames(mat_genus_f))

anno_species <- tax_tbl |>
  mutate(species_grp = map_asv_species[asv_id]) |>
  group_by(species_grp) |>
  summarise(
    domain = safe_mode(domain),
    phylum = safe_mode(phylum),
    class  = safe_mode(class),
    order  = safe_mode(order),
    family = safe_mode(family),
    genus  = safe_mode(genus),
    .groups = "drop"
  ) |>
  mutate(
    feature = species_grp,
    species = species_grp,
    domain = fill_unknown(domain, "domain"),
    phylum = fill_unknown(phylum, "phylum"),
    class  = fill_unknown(class, "class"),
    order  = fill_unknown(order, "order"),
    family = fill_unknown(family, "family"),
    genus  = fill_unknown(genus, "genus")
  ) |>
  select(feature, domain, phylum, class, order, family, genus, species) |>
  filter(feature %in% rownames(mat_species_f))

## ---------- run models ----------
models <- c("unadjusted", "primary", "sensitivity_add_dose")

for (m in models) run_wilcox_suite("ASV",     mat_asv_f,     anno_asv,     meta, m)
for (m in models) run_wilcox_suite("Genus",   mat_genus_f,   anno_genus,   meta, m)
for (m in models) run_wilcox_suite("Species", mat_species_f, anno_species, meta, m)

## ---------- ASV presence/absence supplement ----------
if (RUN_ASV_PRESENCE_SUPPLEMENT) {
  run_asv_presence_supplement(mat_asv, meta, tax_tbl)
}

message("16S DA upgraded done -> ", OUT_ROOT)