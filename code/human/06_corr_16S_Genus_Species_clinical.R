## =========================
## 08_corr_16S_Genus_Species_clinical.R
## 16S clinical association analysis (Genus + Species) with CLR
## - Levels: Genus, Species (separate folders)
## - Covariates: Age, BMI, Disease duration, Constipation, Levodopa daily dose,
##               MDS-UPDRS-III reduction (%)
## - Models: unadjusted / primary / sensitivity_add_dose
##   (via residualization; Spearman for numeric covariates, Wilcoxon for binary)
## - Outputs: full result tables + per-sample reproducible matrices + PDF heatmaps
## Inputs in working dir:
##   metadata.csv
##   table_filtered_w_tax.xls  或  table_f_even_depth_w_tax
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
INPUT_16S <- "table_filtered_w_tax.xls"  #抽平/未抽平

OUT_ROOT <- "results_16S_correlation"  #抽平/未抽平
GROUP_LEVELS <- c("PD_GR", "PD_MR")   # GR reference

# Feature filtering (after aggregation)
MIN_PREVALENCE_ALL <- 0.10
MIN_TOTAL_COUNT <- 10
PSEUDOCOUNT <- 0.5

# Heatmap inclusion: features with min p < 0.05; if none fallback top N
P_CUTOFF <- 0.05
FALLBACK_TOP_N <- 50

# What to run
RUN_WITHIN_GROUP <- TRUE
RUN_INTERACTION  <- TRUE

# Figures
MAKE_PDF <- TRUE
MAKE_PNG <- TRUE
CLUSTER_ROWS <- TRUE
CLUSTER_COLS <- FALSE   # keep sample ordering GR then MR

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
  dir.create(file.path(path, "per_sample"), showWarnings = FALSE)
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

aggregate_counts <- function(mat_feat_by_samp, feature_to_group_named) {
  grp <- feature_to_group_named[rownames(mat_feat_by_samp)]
  rowsum(mat_feat_by_samp, group = grp, reorder = TRUE)
}

filter_features <- function(mat_feat_by_samp) {
  prev <- rowMeans(mat_feat_by_samp > 0)
  tot  <- rowSums(mat_feat_by_samp)
  keep <- (prev >= MIN_PREVALENCE_ALL) & (tot >= MIN_TOTAL_COUNT)
  list(
    mat = mat_feat_by_samp[keep, , drop = FALSE],
    diag = tibble(feature = rownames(mat_feat_by_samp),
                  prev_all = prev,
                  total = tot,
                  keep = keep),
    summary = tibble(n_total=nrow(mat_feat_by_samp), n_keep=sum(keep),
                     min_prev_all=MIN_PREVALENCE_ALL, min_total_count=MIN_TOTAL_COUNT)
  )
}

clr_transform <- function(mat_feat_by_samp, pseudocount = 0.5) {
  x <- mat_feat_by_samp + pseudocount
  x <- sweep(x, 2, colSums(x), "/")
  lx <- log(x)
  sweep(lx, 2, colMeans(lx), "-")
}

residualize_matrix <- function(Y_samp_by_feat, X_design) {
  qrX <- qr(X_design)
  qr.resid(qrX, Y_samp_by_feat)
}

residualize_vector <- function(y, X_design) {
  fit <- lm.fit(X_design, y)
  fit$residuals
}

# Spearman rho + approximate p (t approximation)
spearman_fast <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 6) return(c(rho=NA_real_, p=NA_real_, n=n))
  rho <- suppressWarnings(cor(x, y, method = "spearman"))
  if (is.na(rho)) return(c(rho=NA_real_, p=NA_real_, n=n))
  # t approximation
  rho2 <- min(0.999999, max(-0.999999, rho))
  tval <- rho2 * sqrt((n - 2) / (1 - rho2^2))
  p <- 2 * pt(-abs(tval), df = n - 2)
  c(rho=rho, p=p, n=n)
}

# Wilcoxon + rank-biserial (AUC-based)
wilcox_rb <- function(x, g01) {
  ok <- is.finite(x) & !is.na(g01)
  x <- x[ok]; g01 <- g01[ok]
  if (length(unique(g01)) < 2) return(c(rb=NA_real_, p=NA_real_, n=length(x)))
  x1 <- x[g01 == 1]
  x0 <- x[g01 == 0]
  n1 <- length(x1); n0 <- length(x0)
  if (n1 < 2 || n0 < 2) return(c(rb=NA_real_, p=NA_real_, n=n1+n0))

  wt <- suppressWarnings(wilcox.test(x1, x0, exact = FALSE))
  p <- wt$p.value
  # U via ranks of x1
  z <- c(x1, x0)
  r <- rank(z, ties.method = "average")
  Wrank <- sum(r[seq_len(n1)])
  U <- Wrank - n1*(n1+1)/2
  rb <- (2*U)/(n1*n0) - 1
  rb <- max(-1, min(1, rb))
  c(rb=rb, p=p, n=n1+n0)
}

plot_heatmap_effect <- function(effect_mat, out_pdf, out_png) {
  # effect_mat: features x covariates
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    cr <- CLUSTER_ROWS && nrow(effect_mat) >= 2
    cc <- CLUSTER_COLS && ncol(effect_mat) >= 2
    if (MAKE_PDF) {
      pdf(out_pdf, width = 9, height = max(4, 0.18*nrow(effect_mat)+2))
      pheatmap::pheatmap(effect_mat,
                         cluster_rows = cr,
                         cluster_cols = cc,
                         show_colnames = TRUE,
                         fontsize_row = 8)
      dev.off()
    }
    if (MAKE_PNG) {
      png(out_png, width = 2200, height = max(1400, 40*nrow(effect_mat)+600), res = 300)
      pheatmap::pheatmap(effect_mat,
                         cluster_rows = cr,
                         cluster_cols = cc,
                         show_colnames = TRUE,
                         fontsize_row = 8)
      dev.off()
    }
  } else {
    df <- as.data.frame(effect_mat) |>
      rownames_to_column("feature") |>
      pivot_longer(-feature, names_to="covariate", values_to="effect")
    df$feature <- factor(df$feature, levels = rev(rownames(effect_mat)))
    p <- ggplot(df, aes(covariate, feature, fill = effect)) +
      geom_tile() + theme_bw() +
      labs(x=NULL, y=NULL)
    if (MAKE_PDF) ggsave(out_pdf, p, width=9, height=max(4,0.18*nrow(effect_mat)+2))
    if (MAKE_PNG) ggsave(out_png, p, width=9, height=max(4,0.18*nrow(effect_mat)+2), dpi=300)
  }
}

## ---------- load metadata ----------
meta_raw <- read_csv_smart(META_FILE) |>
  janitor::clean_names(replace = c("\u03bc" = "u"))

sample_col <- find_col(meta_raw, c("^sample$"))
group_col  <- find_col(meta_raw, c("^group$"))

age_col  <- find_col(meta_raw, c("^age$"))
bmi_col  <- find_col(meta_raw, c("^bmi"))
dur_col  <- find_col(meta_raw, c("disease_duration"))
con_col  <- find_col(meta_raw, c("constipation"))
dose_col <- find_col(meta_raw, c("levodopa.*dose|levodopa_daily_dose"))

# NEW covariate: MDS-UPDRS-III reduction (%)
updrs_red_col <- find_col(meta_raw, c("mds.*updrs.*iii.*reduction", "updrs.*iii.*reduction", "reduction.*updrs.*iii"))

stopifnot(!is.na(sample_col), !is.na(group_col))
if (is.na(updrs_red_col)) {
  stop("Cannot find column for MDS-UPDRS-III reduction (%). Please check metadata column names after clean_names().")
}

parse_percent <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "%", "")
  suppressWarnings(as.numeric(x))
}

meta <- meta_raw |>
  mutate(sample_id = normalize_id(.data[[sample_col]]),
         group = .data[[group_col]]) |>
  filter(!is.na(sample_id) & sample_id != "") |>
  filter(group %in% GROUP_LEVELS) |>
  transmute(
    sample_id,
    group = factor(group, levels = GROUP_LEVELS),
    age = suppressWarnings(as.numeric(.data[[age_col]])),
    bmi = suppressWarnings(as.numeric(.data[[bmi_col]])),
    disease_duration = suppressWarnings(as.numeric(.data[[dur_col]])),
    constipation = factor(as.integer(.data[[con_col]]), levels = c(0,1), labels = c("No","Yes")),
    levodopa_dose = suppressWarnings(as.numeric(.data[[dose_col]])),
    updrs3_reduction_pct = parse_percent(.data[[updrs_red_col]])
  )

stopifnot(nrow(meta) == 40)
meta <- meta |> arrange(group, sample_id)
sample_order <- meta$sample_id

dir.create(OUT_ROOT, showWarnings = FALSE)
write_tsv0(meta, file.path(OUT_ROOT, "00_manifest_metadata_model_40.tsv"))
write_tsv0(tibble(
  age_col=age_col, bmi_col=bmi_col, dur_col=dur_col, con_col=con_col, dose_col=dose_col, updrs_red_col=updrs_red_col
), file.path(OUT_ROOT, "00_manifest_metadata_column_mapping.tsv"))

## ---------- read 16S table (often TSV with .xls extension) ----------
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
  mutate(
    genus = fill_unknown(genus, "genus"),
    species = fill_unknown(species, "species"),
    family = fill_unknown(family, "family"),
    phylum = fill_unknown(phylum, "phylum")
  )

write_tsv0(tax_tbl, file.path(OUT_ROOT, "00_manifest_ASV_taxonomy_parsed.tsv"))

## ---------- build Genus & Species matrices ----------
map_asv_genus <- setNames(tax_tbl$genus, tax_tbl$asv_id)
mat_genus <- aggregate_counts(mat_asv, map_asv_genus)

map_asv_species <- setNames(tax_tbl$species, tax_tbl$asv_id)
mat_species <- aggregate_counts(mat_asv, map_asv_species)

# filter
fg <- filter_features(mat_genus)
fs <- filter_features(mat_species)
mat_genus_f <- fg$mat
mat_species_f <- fs$mat

write_tsv0(fg$summary, file.path(OUT_ROOT, "00_manifest_filter_summary_Genus.tsv"))
write_tsv0(fs$summary, file.path(OUT_ROOT, "00_manifest_filter_summary_Species.tsv"))
write_tsv0(fg$diag, file.path(OUT_ROOT, "00_manifest_filter_diagnostics_Genus.tsv"))
write_tsv0(fs$diag, file.path(OUT_ROOT, "00_manifest_filter_diagnostics_Species.tsv"))

# annotations for genus/species
anno_genus <- tax_tbl |>
  mutate(genus_grp = map_asv_genus[asv_id]) |>
  group_by(genus_grp) |>
  summarise(
    phylum = safe_mode(phylum),
    family = safe_mode(family),
    .groups="drop"
  ) |>
  mutate(feature = genus_grp) |>
  select(feature, phylum, family) |>
  filter(feature %in% rownames(mat_genus_f))

anno_species <- tax_tbl |>
  mutate(species_grp = map_asv_species[asv_id]) |>
  group_by(species_grp) |>
  summarise(
    genus = safe_mode(genus),
    phylum = safe_mode(phylum),
    family = safe_mode(family),
    .groups="drop"
  ) |>
  mutate(feature = species_grp) |>
  select(feature, genus, phylum, family) |>
  filter(feature %in% rownames(mat_species_f))

## ---------- correlation runner ----------
COVARS <- tibble::tribble(
  ~covar, ~type,
  "age", "numeric",
  "bmi", "numeric",
  "disease_duration", "numeric",
  "constipation", "binary",
  "levodopa_dose", "numeric",
  "updrs3_reduction_pct", "numeric"
)

build_adjust_X <- function(meta, model_name, tested_covar) {
  # returns design matrix for adjustment (excluding tested_covar)
  # primary adjustment set: group + age + bmi + disease_duration + constipation
  # sensitivity adds levodopa_dose
  base <- c("group", "age", "bmi", "disease_duration", "constipation")
  if (model_name == "sensitivity_add_dose") base <- c(base, "levodopa_dose")

  adj <- setdiff(base, tested_covar)
  if (model_name == "unadjusted") return(NULL)

  # if nothing left beyond intercept, return intercept-only matrix
  if (length(adj) == 0) return(model.matrix(~ 1, data = meta))
  f <- as.formula(paste("~", paste(adj, collapse = " + ")))
  model.matrix(f, data = meta)
}

run_one_level <- function(level_name, mat_counts, feat_anno) {
  level_root <- file.path(OUT_ROOT, level_name)
  dir_init(level_root)

  # order samples
  mat_counts <- mat_counts[, sample_order, drop=FALSE]
  # CLR (features x samples)
  clr_feat_by_samp <- clr_transform(mat_counts, pseudocount = PSEUDOCOUNT)

  # save base per-sample matrices
  write_tsv0(as.data.frame(clr_feat_by_samp) |> rownames_to_column("feature"),
             file.path(level_root, "per_sample", "CLR_features_by_sample.tsv"))
  write_tsv0(as.data.frame(mat_counts) |> rownames_to_column("feature"),
             file.path(level_root, "per_sample", "Counts_features_by_sample.tsv"))
  write_tsv0(feat_anno, file.path(level_root, "00_manifest", "feature_annotation.tsv"))
  write_tsv0(COVARS, file.path(level_root, "00_manifest", "covariates.tsv"))

  models <- c("unadjusted","primary","sensitivity_add_dose")

  for (mname in models) {
    out_dir <- file.path(level_root, mname)
    dir_init(out_dir)
    write_tsv0(tibble(model=mname), file.path(out_dir, "00_manifest", "model.tsv"))

    # working matrix: samples x features
    Y0 <- t(clr_feat_by_samp)
    rownames(Y0) <- colnames(clr_feat_by_samp)
    colnames(Y0) <- rownames(clr_feat_by_samp)

    all_res <- list()

    n_feat <- ncol(Y0)  # number of features
    eff_mat <- matrix(NA_real_, nrow = n_feat, ncol = nrow(COVARS))
    rownames(eff_mat) <- colnames(Y0)   # feature names
    colnames(eff_mat) <- COVARS$covar

    p_mat <- eff_mat
    q_mat <- eff_mat

    for (i in seq_len(nrow(COVARS))) {
      cv <- COVARS$covar[i]
      tp <- COVARS$type[i]

      Xadj <- build_adjust_X(meta, mname, cv)

      # feature residuals (if adjusted)
      Y <- Y0
      if (!is.null(Xadj) && mname != "unadjusted") {
        Y <- residualize_matrix(Y0, Xadj)
      }

      # numeric covariate residuals (for partial spearman)
      cov_raw <- meta[[cv]]
      cov_use <- cov_raw

      if (tp == "numeric" && !is.null(Xadj) && mname != "unadjusted") {
        cov_use <- residualize_vector(cov_raw, Xadj)
      }

      if (tp == "numeric") {
        # compute per feature
        for (f in colnames(Y)) {
          ss <- spearman_fast(Y[, f], cov_use)
          all_res[[length(all_res)+1]] <- tibble(
            feature = f, covariate = cv, test = "spearman",
            effect = ss["rho"], p_value = ss["p"], n = ss["n"]
          )
          eff_mat[f, cv] <- as.numeric(ss["rho"])
          p_mat[f, cv] <- as.numeric(ss["p"])
        }
      } else {
        # binary (constipation): compare residualized feature by group (0/1)
        g01 <- ifelse(meta[[cv]] == "Yes", 1L, 0L)
        for (f in colnames(Y)) {
          ww <- wilcox_rb(Y[, f], g01)
          all_res[[length(all_res)+1]] <- tibble(
            feature = f, covariate = cv, test = "wilcox",
            effect = ww["rb"], p_value = ww["p"], n = ww["n"]
          )
          eff_mat[f, cv] <- as.numeric(ww["rb"])
          p_mat[f, cv] <- as.numeric(ww["p"])
        }
      }

      # reproducible per-sample data used for this covariate/model
      # feature values used:
      write_tsv0(as.data.frame(t(Y)) |> rownames_to_column("feature"),
                 file.path(out_dir, "per_sample", paste0("feature_values_used_for_", cv, ".tsv")))
      # covariate used (numeric residual or raw)
      if (tp == "numeric") {
        write_tsv0(tibble(sample_id = meta$sample_id, covariate = cv, value_used = cov_use),
                   file.path(out_dir, "per_sample", paste0("covariate_values_used_for_", cv, ".tsv")))
      }
    }

    res_df <- bind_rows(all_res) |>
      group_by(covariate) |>
      mutate(q_value = p.adjust(p_value, method="BH")) |>
      ungroup() |>
      left_join(feat_anno, by="feature")

    write_tsv0(res_df, file.path(out_dir, "tables", "association_overall_long.tsv"))

    # matrices (effects/p/q)
    # q matrix from res_df:
    q_tmp <- res_df |> select(feature, covariate, q_value) |>
      pivot_wider(names_from=covariate, values_from=q_value)
    q_mat2 <- as.data.frame(q_tmp); rownames(q_mat2) <- q_mat2$feature; q_mat2$feature <- NULL
    q_mat[rownames(q_mat2), colnames(q_mat2)] <- as.matrix(q_mat2)

    write_tsv0(as.data.frame(eff_mat) |> rownames_to_column("feature"),
               file.path(out_dir, "tables", "association_effect_matrix.tsv"))
    write_tsv0(as.data.frame(p_mat) |> rownames_to_column("feature"),
               file.path(out_dir, "tables", "association_p_matrix.tsv"))
    write_tsv0(as.data.frame(q_mat) |> rownames_to_column("feature"),
               file.path(out_dir, "tables", "association_q_matrix.tsv"))

    # heatmap selection
    min_p <- apply(p_mat, 1, function(x) suppressWarnings(min(x, na.rm=TRUE)))
    sel <- names(min_p)[is.finite(min_p) & min_p < P_CUTOFF]
    if (length(sel) == 0) {
      sel <- names(sort(min_p))[seq_len(min(FALLBACK_TOP_N, length(min_p)))]
      writeLines(paste0("No features with min p<", P_CUTOFF, "; fallback top ", length(sel), " by min p."),
                 con = file.path(out_dir, "logs", "heatmap_selection_note.txt"))
    } else {
      sel <- sel[order(min_p[sel])]
      if (length(sel) > FALLBACK_TOP_N) sel <- sel[seq_len(FALLBACK_TOP_N)]
      writeLines(paste0("Selected features: min p<", P_CUTOFF, " (n=", length(sel), "), capped at ", length(sel), "."),
                 con = file.path(out_dir, "logs", "heatmap_selection_note.txt"))
    }
    write_tsv0(tibble(feature = sel, min_p = min_p[sel]),
               file.path(out_dir, "tables", "heatmap_selected_features.tsv"))

    eff_sel <- eff_mat[sel, , drop=FALSE]
    # replace NA with 0 for plotting aesthetics? keep NA
    if (MAKE_PDF || MAKE_PNG) {
      plot_heatmap_effect(eff_sel,
                          out_pdf = file.path(out_dir, "figures", paste0("heatmap_effect_", level_name, "_", mname, ".pdf")),
                          out_png = file.path(out_dir, "figures", paste0("heatmap_effect_", level_name, "_", mname, ".png")))
    }

    # within-group (optional): unadjusted only by default; here we do same model scheme but without group adjustment
    if (RUN_WITHIN_GROUP) {
      for (gname in GROUP_LEVELS) {
        idx <- which(meta$group == gname)
        meta_g <- meta[idx, , drop=FALSE]
        Yg0 <- t(clr_feat_by_samp[, idx, drop=FALSE])
        rownames(Yg0) <- meta_g$sample_id
        colnames(Yg0) <- rownames(clr_feat_by_samp)

        allg <- list()
        for (i in seq_len(nrow(COVARS))) {
          cv <- COVARS$covar[i]
          tp <- COVARS$type[i]

          # adjustment set within-group (no group term)
          base <- c("age","bmi","disease_duration","constipation")
          if (mname == "sensitivity_add_dose") base <- c(base, "levodopa_dose")
          adj <- setdiff(base, cv)
          Xadj <- if (mname == "unadjusted") NULL else {
            if (length(adj) == 0) model.matrix(~1, data = meta_g)
            else model.matrix(as.formula(paste("~", paste(adj, collapse=" + "))), data = meta_g)
          }

          Yg <- Yg0
          if (!is.null(Xadj) && mname != "unadjusted") Yg <- residualize_matrix(Yg0, Xadj)

          cov_raw <- meta_g[[cv]]
          cov_use <- cov_raw
          if (tp == "numeric" && !is.null(Xadj) && mname != "unadjusted") cov_use <- residualize_vector(cov_raw, Xadj)

          if (tp == "numeric") {
            for (f in colnames(Yg)) {
              ss <- spearman_fast(Yg[, f], cov_use)
              allg[[length(allg)+1]] <- tibble(
                group = gname, model = mname, feature = f, covariate = cv,
                test = "spearman", effect = ss["rho"], p_value = ss["p"], n = ss["n"]
              )
            }
          } else {
            g01 <- ifelse(meta_g[[cv]] == "Yes", 1L, 0L)
            for (f in colnames(Yg)) {
              ww <- wilcox_rb(Yg[, f], g01)
              allg[[length(allg)+1]] <- tibble(
                group = gname, model = mname, feature = f, covariate = cv,
                test = "wilcox", effect = ww["rb"], p_value = ww["p"], n = ww["n"]
              )
            }
          }
        }
        gdf <- bind_rows(allg) |>
          group_by(covariate) |>
          mutate(q_value = p.adjust(p_value, method="BH")) |>
          ungroup() |>
          left_join(feat_anno, by="feature")
        write_tsv0(gdf, file.path(out_dir, "tables", paste0("association_within_", gname, "_long.tsv")))
      }
    }

    # interaction (optional): linear model on CLR scale
    if (RUN_INTERACTION) {
      # use CLR (not residual) for interaction models
      Yclr <- t(clr_feat_by_samp)
      idx_gr <- which(meta$group == "PD_GR")
      idx_mr <- which(meta$group == "PD_MR")

      inter_all <- list()
      for (i in seq_len(nrow(COVARS))) {
        cv <- COVARS$covar[i]
        # build adjustment set per model (same as earlier) but include main effects
        base <- c("age","bmi","disease_duration","constipation")
        if (mname == "sensitivity_add_dose") base <- c(base, "levodopa_dose")
        others <- setdiff(base, cv)

        rhs <- paste0("group * ", cv)
        if (length(others) > 0) rhs <- paste(rhs, paste(others, collapse=" + "), sep=" + ")
        fml <- as.formula(paste0("y ~ ", rhs))

        for (f in colnames(Yclr)) {
          df <- meta |> mutate(y = Yclr[, f])
          fit <- try(lm(fml, data = df), silent = TRUE)
          if (inherits(fit, "try-error")) next
          coefs <- summary(fit)$coefficients
          # interaction term name patterns
          iterm <- grep(paste0("^groupPD_MR:", cv, "$|^", cv, ":groupPD_MR$"), rownames(coefs), value=TRUE)
          if (length(iterm) == 0) next
          it <- iterm[1]
          inter_all[[length(inter_all)+1]] <- tibble(
            feature = f, covariate = cv,
            beta_interaction = coefs[it, "Estimate"],
            p_value = coefs[it, "Pr(>|t|)"]
          )
        }
      }
      inter_df <- bind_rows(inter_all) |>
        group_by(covariate) |>
        mutate(q_value = p.adjust(p_value, method="BH")) |>
        ungroup() |>
        left_join(feat_anno, by="feature")
      write_tsv0(inter_df, file.path(out_dir, "tables", "association_interaction_long.tsv"))
    }

    message("Done correlation: ", level_name, " | ", mname)
  }
}

## ---------- run Genus + Species (separate folders) ----------
run_one_level("Genus", mat_genus_f, anno_genus)
run_one_level("Species", mat_species_f, anno_species)

message("All 16S correlation analyses done -> ", OUT_ROOT)