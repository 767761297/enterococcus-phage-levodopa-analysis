#!/usr/bin/env Rscript
## =========================================================
## alpha-metrics vs clinical correlation + UPDRS scatter
## - Heatmap: X=clinical variables, Y=alpha metrics
## - Scatter: for ALL alpha metrics vs UPDRS (overall + groupwise)
## - Models for correlation tables: unadjusted / primary / sensitivity_add_dose
##   (same logic as 06_corr_16S_Genus_Species_clinical.R)
## - Residualization for scatter: same as 06_scatter_16S_updrs3_Enterococcus.R
##   (NO group term; within-group residualization for groupwise)
## =========================================================

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
  library(pheatmap)
})

## ------------------------
## 0) USER EDITABLE PARAMS
## ------------------------
# Input tables (per-sample alpha with metadata merged preferred)
# 16S Enterococcus within-genus alpha (species level; includes Chao1)
INPUT_16S_ALPHA  <- "alpha_per_sample_Enterococcus_species.tsv"
# virome Enterococcus-host Host_species alpha (TPM-based; no Chao1)
INPUT_VIROME_ALPHA <- "alpha_per_sample_EnterococcusHostSpeciesPhage.tsv"

# If alpha tables don't include clinical columns, provide metadata.csv for merging
META_FILE <- "metadata.csv"

# Output root
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_ROOT <- paste0("results_alpha_vs_clinical_", timestamp)

# Groups / colors consistent with your old scripts
GROUP_LEVELS <- c("PD_GR","PD_MR")
COL_MAP <- c(PD_GR = "#7CAE00", PD_MR = "#F8766D")
LAB_MAP <- c(PD_GR = "GR", PD_MR = "MR")

MODELS <- c("unadjusted","primary","sensitivity_add_dose")

MAKE_PDF <- TRUE
MAKE_PNG <- TRUE

# Heatmap significance symbols based on q-value
Q_CUTS <- c(0.05, 0.10)  # marks: ** for <0.05, * for <0.10

## ---------------------------------
## 1) Helpers (borrowed from old style)
## ---------------------------------
dir_init <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(path, "00_manifest"), showWarnings = FALSE)
  dir.create(file.path(path, "tables"), showWarnings = FALSE)
  dir.create(file.path(path, "figures_pdf"), showWarnings = FALSE)
  dir.create(file.path(path, "figures_png"), showWarnings = FALSE)
  dir.create(file.path(path, "per_sample"), showWarnings = FALSE)
  dir.create(file.path(path, "logs"), showWarnings = FALSE)
}

write_tsv0 <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(x, path)
}

read_tsv_smart <- function(path) {
  # TSV/CSV with encoding fallback
  tryCatch(
    readr::read_tsv(path, show_col_types = FALSE, locale = locale(encoding = "UTF-8")),
    error = function(e) readr::read_tsv(path, show_col_types = FALSE, locale = locale(encoding = "GB18030"))
  )
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

parse_percent <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "%", "")
  suppressWarnings(as.numeric(x))
}

# residualize vector on design matrix
residualize_vector <- function(y, X) {
  ok <- complete.cases(X) & is.finite(y)
  out <- rep(NA_real_, length(y))
  if (sum(ok) < 6) return(out)
  fit <- lm.fit(X[ok, , drop=FALSE], y[ok])
  out[ok] <- fit$residuals
  out
}
# residualize matrix (samples x features)
residualize_matrix <- function(Y, X) {
  # Y: samples x features matrix (may contain NA/NaN/Inf)
  # X: design matrix (samples x p)
  stopifnot(nrow(Y) == nrow(X))

  okX <- complete.cases(X)
  out <- matrix(NA_real_, nrow=nrow(Y), ncol=ncol(Y), dimnames=dimnames(Y))
  if (sum(okX) < 6) return(out)

  for (j in seq_len(ncol(Y))) {
    y <- Y[, j]
    ok <- okX & is.finite(y)
    if (sum(ok) < 6) next
    fit <- lm.fit(X[ok, , drop=FALSE], y[ok])
    out[ok, j] <- fit$residuals
  }
  out
}
# Spearman rho + t approximation p (same style as your old script)
spearman_fast <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 6) return(c(rho=NA_real_, p=NA_real_, n=n))
  rho <- suppressWarnings(cor(x, y, method="spearman"))
  if (is.na(rho)) return(c(rho=NA_real_, p=NA_real_, n=n))
  rho2 <- min(0.999999, max(-0.999999, rho))
  tval <- rho2 * sqrt((n - 2) / (1 - rho2^2))
  p <- 2 * pt(-abs(tval), df = n - 2)
  c(rho=rho, p=p, n=n)
}

# Wilcoxon + rank-biserial (rb) effect
# g01: 0/1; interpret effect as "1 vs 0" (positive => higher in group=1)
wilcox_rb <- function(x, g01) {
  ok <- is.finite(x) & !is.na(g01)
  x <- x[ok]; g01 <- g01[ok]
  n <- length(x)
  if (n < 6) return(c(rb=NA_real_, p=NA_real_, n=n))
  g01 <- as.integer(g01)
  if (length(unique(g01)) < 2) return(c(rb=NA_real_, p=NA_real_, n=n))
  # ensure group1=1, group0=0
  x1 <- x[g01==1]; x0 <- x[g01==0]
  n1 <- length(x1); n0 <- length(x0)
  if (n1 < 2 || n0 < 2) return(c(rb=NA_real_, p=NA_real_, n=n))
  wt <- suppressWarnings(wilcox.test(x1, x0, exact=FALSE))
  # wilcox.test returns W = sum of ranks for x1
  W <- as.numeric(wt$statistic)
  U <- W - n1*(n1+1)/2
  AUC <- U / (n1*n0)
  rb <- 2*AUC - 1
  c(rb=rb, p=wt$p.value, n=n1+n0)
}

## -------------------------------
## 2) Load metadata + normalize columns
## -------------------------------
meta_raw <- read_csv_smart(META_FILE) |> janitor::clean_names(replace = c("\u03bc"="u"))

sample_col <- find_col(meta_raw, c("^sample$"))
group_col  <- find_col(meta_raw, c("^group$"))
age_col    <- find_col(meta_raw, c("^age$"))
bmi_col    <- find_col(meta_raw, c("^bmi"))
dur_col    <- find_col(meta_raw, c("disease_duration"))
con_col    <- find_col(meta_raw, c("constipation"))
dose_col   <- find_col(meta_raw, c("levodopa.*dose|levodopa_daily_dose"))
updrs_col  <- find_col(meta_raw, c("mds.*updrs.*iii.*reduction", "updrs.*iii.*reduction", "reduction.*updrs.*iii"))

stopifnot(!is.na(sample_col), !is.na(group_col), !is.na(updrs_col))

meta <- meta_raw |>
  mutate(
    sample_id = normalize_id(.data[[sample_col]]),
    group = .data[[group_col]]
  ) |>
  filter(group %in% GROUP_LEVELS) |>
  transmute(
    sample_id,
    group = factor(group, levels = GROUP_LEVELS),
    age = suppressWarnings(as.numeric(.data[[age_col]])),
    bmi = suppressWarnings(as.numeric(.data[[bmi_col]])),
    disease_duration = suppressWarnings(as.numeric(.data[[dur_col]])),
    constipation = factor(as.integer(.data[[con_col]]), levels=c(0,1), labels=c("No","Yes")),
    levodopa_dose = suppressWarnings(as.numeric(.data[[dose_col]])),
    updrs3_reduction_pct = parse_percent(.data[[updrs_col]])
  ) |>
  arrange(group, sample_id)

# sample order used everywhere
sample_order <- meta$sample_id

## covariates definition (X-axis for heatmap)
COVARS <- tibble::tribble(
  ~covar, ~type,
  "age", "numeric",
  "bmi", "numeric",
  "disease_duration", "numeric",
  "constipation", "binary",
  "levodopa_dose", "numeric",
  "updrs3_reduction_pct", "numeric"
)

## -------------------------------
## 3) Load alpha tables and standardize
## -------------------------------
load_alpha <- function(path, dataset_name) {
  df <- read_tsv_smart(path) |> janitor::clean_names()
  # find sample id col
  sid <- find_col(df, c("^sample_id$", "^sample$"))
  stopifnot(!is.na(sid))
  df <- df |> mutate(sample_id = normalize_id(.data[[sid]]))
  # merge metadata if missing key clinical fields
  need_merge <- any(!c("group","age","bmi","disease_duration","constipation","levodopa_dose","updrs3_reduction_pct") %in% names(df))
  if (need_merge) {
    df <- df |> select(-matches("^group$|^age$|^bmi$|^disease_duration$|^constipation$|^levodopa_dose$|^updrs3_reduction_pct$"))
    df <- df |> left_join(meta, by="sample_id")
  } else {
    # ensure group factor levels & types
    df <- df |>
      mutate(
        group = factor(as.character(group), levels = GROUP_LEVELS),
        constipation = factor(as.character(constipation), levels=c("No","Yes")),
        updrs3_reduction_pct = suppressWarnings(as.numeric(updrs3_reduction_pct)),
        age = suppressWarnings(as.numeric(age)),
        bmi = suppressWarnings(as.numeric(bmi)),
        disease_duration = suppressWarnings(as.numeric(disease_duration)),
        levodopa_dose = suppressWarnings(as.numeric(levodopa_dose))
      )
  }
  df <- df |> filter(sample_id %in% sample_order) |> arrange(group, sample_id)
  df$dataset <- dataset_name
  df
}

alpha16 <- load_alpha(INPUT_16S_ALPHA, "16S_Enterococcus_within_genus")
alphav <- load_alpha(INPUT_VIROME_ALPHA, "Virome_EnterococcusHostSpecies")

## alpha metric lists (Y-axis for heatmap)
METRICS_16S <- c("observed","chao1","shannon","exp_shannon","simpson","inv_simpson","pielou_evenness")
METRICS_VIR <- c("observed_gt0","observed_gt1","shannon","exp_shannon","simpson","inv_simpson","pielou_evenness")

## -------------------------------
## 4) Build adjustment design matrix for correlation tables
##    (same as 06_corr_16S_Genus_Species_clinical.R)
## -------------------------------
build_adjust_X <- function(meta, model_name, tested_covar) {
  base <- c("group","age","bmi","disease_duration","constipation")
  if (model_name == "sensitivity_add_dose") base <- c(base, "levodopa_dose")
  adj <- setdiff(base, tested_covar)
  if (model_name == "unadjusted") return(NULL)
  if (length(adj) == 0) return(model.matrix(~ 1, data = meta))
  f <- as.formula(paste("~", paste(adj, collapse = " + ")))
  model.matrix(f, data = meta)
}

## -------------------------------
## 5) Correlation runner (overall + within-group + interaction) for alpha metrics
## -------------------------------
run_alpha_clinical <- function(df_alpha, metrics, dataset_tag, out_root) {

  dir_init(out_root)
  write_tsv0(meta, file.path(out_root, "00_manifest", "metadata_model_40.tsv"))
  write_tsv0(COVARS, file.path(out_root, "00_manifest", "covariates.tsv"))
  write_tsv0(tibble(metric=metrics), file.path(out_root, "00_manifest", "alpha_metrics.tsv"))

  # base feature matrix: samples x features
  # use sample_order to align
  Y0 <- df_alpha |>
    select(sample_id, all_of(metrics)) |>
    filter(sample_id %in% sample_order) |>
    arrange(match(sample_id, sample_order)) |>
    column_to_rownames("sample_id") |>
    as.matrix()

Y0[!is.finite(Y0)] <- NA_real_

  # For interaction: use raw alpha (not residual) like old script "use CLR not residual"
  # Here we use raw alpha metrics.
  for (mname in MODELS) {
    out_dir <- file.path(out_root, mname)
    dir_init(out_dir)
    write_tsv0(tibble(model=mname), file.path(out_dir, "00_manifest", "model.tsv"))

    eff_mat <- matrix(NA_real_, nrow=ncol(Y0), ncol=nrow(COVARS),
                      dimnames=list(colnames(Y0), COVARS$covar))
    p_mat <- eff_mat
    q_mat <- eff_mat
    all_res <- list()

    # loop covariates
    for (i in seq_len(nrow(COVARS))) {
      cv <- COVARS$covar[i]
      tp <- COVARS$type[i]
      Xadj <- build_adjust_X(meta, mname, cv)

      # feature residuals
      Y <- Y0
      if (!is.null(Xadj) && mname != "unadjusted") {
        Y <- residualize_matrix(Y0, Xadj)
      }

      cov_raw <- meta[[cv]]
      cov_use <- cov_raw

      if (tp == "numeric" && !is.null(Xadj) && mname != "unadjusted") {
        cov_use <- residualize_vector(cov_raw, Xadj)
      }

      # save per-sample used values (for reproducibility)
      write_tsv0(
        as.data.frame(Y) |> rownames_to_column("sample_id"),
        file.path(out_dir, "per_sample", paste0("feature_values_used_for_", cv, ".tsv"))
      )
      write_tsv0(
        tibble(sample_id = meta$sample_id, covariate=cv, value_used = cov_use),
        file.path(out_dir, "per_sample", paste0("covariate_values_used_for_", cv, ".tsv"))
      )

      if (tp == "numeric") {
        for (f in colnames(Y)) {
          ss <- spearman_fast(Y[, f], cov_use)
          all_res[[length(all_res)+1]] <- tibble(
            dataset = dataset_tag, model = mname,
            feature = f, covariate = cv, test = "spearman",
            effect = as.numeric(ss["rho"]),
            p_value = as.numeric(ss["p"]),
            n = as.integer(ss["n"])
          )
          eff_mat[f, cv] <- as.numeric(ss["rho"])
          p_mat[f, cv] <- as.numeric(ss["p"])
        }
      } else {
        # binary: constipation
        g01 <- ifelse(meta[[cv]] == "Yes", 1L, 0L)
        for (f in colnames(Y)) {
          ww <- wilcox_rb(Y[, f], g01)
          all_res[[length(all_res)+1]] <- tibble(
            dataset = dataset_tag, model = mname,
            feature = f, covariate = cv, test = "wilcox",
            effect = as.numeric(ww["rb"]),
            p_value = as.numeric(ww["p"]),
            n = as.integer(ww["n"])
          )
          eff_mat[f, cv] <- as.numeric(ww["rb"])
          p_mat[f, cv] <- as.numeric(ww["p"])
        }
      }
    }

    res_df <- bind_rows(all_res) |>
      group_by(covariate) |>
      mutate(q_value = p.adjust(p_value, method="BH")) |>
      ungroup()

    # build q_mat
    q_wide <- res_df |>
      select(feature, covariate, q_value) |>
      pivot_wider(names_from=covariate, values_from=q_value)
    q_mat[rownames(q_mat), colnames(q_mat)] <- as.matrix(q_wide |> column_to_rownames("feature"))[rownames(q_mat), colnames(q_mat)]

    write_tsv0(res_df, file.path(out_dir, "tables", "association_overall_long.tsv"))
    write_tsv0(as.data.frame(eff_mat) |> rownames_to_column("alpha_metric"),
               file.path(out_dir, "tables", "association_effect_matrix.tsv"))
    write_tsv0(as.data.frame(p_mat) |> rownames_to_column("alpha_metric"),
               file.path(out_dir, "tables", "association_p_matrix.tsv"))
    write_tsv0(as.data.frame(q_mat) |> rownames_to_column("alpha_metric"),
               file.path(out_dir, "tables", "association_q_matrix.tsv"))

    ## Heatmap (X=clinical variables, Y=alpha metrics)
    # significance symbols based on q
    sig <- matrix("", nrow=nrow(q_mat), ncol=ncol(q_mat),
                  dimnames=dimnames(q_mat))
    sig[q_mat < Q_CUTS[2]] <- "*"
    sig[q_mat < Q_CUTS[1]] <- "**"

    hm_title <- paste0(dataset_tag, " | ", mname, " | effect (rho/rb); X=clinical, Y=alpha")

    if (MAKE_PDF) {
      pdf(file.path(out_dir, "figures_pdf", "heatmap_alpha_vs_clinical_effect.pdf"),
          width = 8.5, height = max(3.5, 0.3*nrow(eff_mat) + 2.5))
      pheatmap(
        eff_mat,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        display_numbers = sig,
        number_color = "black",
        main = hm_title,
        fontsize = 11
      )
      dev.off()
    }
    if (MAKE_PNG) {
      png(file.path(out_dir, "figures_png", "heatmap_alpha_vs_clinical_effect.png"),
          width = 1300, height = max(700, 130*nrow(eff_mat)), res = 150)
      pheatmap(
        eff_mat,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        display_numbers = sig,
        number_color = "black",
        main = hm_title,
        fontsize = 11
      )
      dev.off()
    }

    ## Within-group correlations (tables only, like old script)
    for (gname in GROUP_LEVELS) {
      meta_g <- meta |> filter(group == gname)
      if (nrow(meta_g) < 6) next
      Y0g <- Y0[meta_g$sample_id, , drop=FALSE]

      allg <- list()
      for (i in seq_len(nrow(COVARS))) {
        cv <- COVARS$covar[i]
        tp <- COVARS$type[i]

        # within group adjustment base excludes group
        base <- c("age","bmi","disease_duration","constipation")
        if (mname == "sensitivity_add_dose") base <- c(base, "levodopa_dose")
        adj <- setdiff(base, cv)
        Xadjg <- if (mname == "unadjusted") NULL else {
          if (length(adj)==0) model.matrix(~1, data=meta_g) else model.matrix(as.formula(paste("~", paste(adj, collapse=" + "))), data=meta_g)
        }

        Yg <- Y0g
        if (!is.null(Xadjg) && mname != "unadjusted") Yg <- residualize_matrix(Y0g, Xadjg)

        cov_raw <- meta_g[[cv]]
        cov_use <- cov_raw
        if (tp=="numeric" && !is.null(Xadjg) && mname != "unadjusted") cov_use <- residualize_vector(cov_raw, Xadjg)

        if (tp=="numeric") {
          for (f in colnames(Yg)) {
            ss <- spearman_fast(Yg[, f], cov_use)
            allg[[length(allg)+1]] <- tibble(
              group = gname, dataset=dataset_tag, model=mname,
              feature=f, covariate=cv, test="spearman",
              effect=as.numeric(ss["rho"]), p_value=as.numeric(ss["p"]), n=as.integer(ss["n"])
            )
          }
        } else {
          g01 <- ifelse(meta_g[[cv]] == "Yes", 1L, 0L)
          for (f in colnames(Yg)) {
            ww <- wilcox_rb(Yg[, f], g01)
            allg[[length(allg)+1]] <- tibble(
              group = gname, dataset=dataset_tag, model=mname,
              feature=f, covariate=cv, test="wilcox",
              effect=as.numeric(ww["rb"]), p_value=as.numeric(ww["p"]), n=as.integer(ww["n"])
            )
          }
        }
      }
      gdf <- bind_rows(allg) |>
        group_by(covariate) |>
        mutate(q_value = p.adjust(p_value, method="BH")) |>
        ungroup()
      write_tsv0(gdf, file.path(out_dir, "tables", paste0("association_within_", gname, "_long.tsv")))
    }

    ## Interaction (numeric covariates only; like old script; tables only)
    inter_all <- list()
    num_covs <- COVARS |> filter(type=="numeric") |> pull(covar)
    for (cv in num_covs) {
      base <- c("age","bmi","disease_duration","constipation")
      if (mname == "sensitivity_add_dose") base <- c(base, "levodopa_dose")
      others <- setdiff(base, cv)

      rhs <- paste0("group * ", cv)
      if (length(others) > 0) rhs <- paste(rhs, paste(others, collapse=" + "), sep=" + ")
      fml <- as.formula(paste0("y ~ ", rhs))

      for (f in colnames(Y0)) {
        dfm <- meta |> mutate(y = Y0[, f])
        fit <- try(lm(fml, data=dfm), silent=TRUE)
        if (inherits(fit, "try-error")) next
        coefs <- summary(fit)$coefficients
        iterm <- grep(paste0("^groupPD_MR:", cv, "$|^", cv, ":groupPD_MR$"), rownames(coefs), value=TRUE)
        if (length(iterm)==0) next
        it <- iterm[1]
        inter_all[[length(inter_all)+1]] <- tibble(
          dataset=dataset_tag, model=mname,
          feature=f, covariate=cv,
          beta_interaction = coefs[it, "Estimate"],
          p_value = coefs[it, "Pr(>|t|)"]
        )
      }
    }
    if (length(inter_all) > 0) {
      inter_df <- bind_rows(inter_all) |>
        group_by(covariate) |>
        mutate(q_value = p.adjust(p_value, method="BH")) |>
        ungroup()
      write_tsv0(inter_df, file.path(out_dir, "tables", "association_interaction_long.tsv"))
    }

    message("Done correlations: ", dataset_tag, " | ", mname)
  }
}

## -------------------------------
## 6) Scatter plots: ALL metrics vs UPDRS (overall + groupwise)
##    (same residualization logic as 06_scatter_16S_updrs3_Enterococcus.R)
## -------------------------------
# residualize vector on covariates (NO group term)
resid_on_covars <- function(df, value, covars) {
  keep <- complete.cases(df[, covars, drop=FALSE]) & is.finite(value)
  out <- rep(NA_real_, length(value))
  if (sum(keep) < 6) return(out)
  X <- model.matrix(as.formula(paste("~", paste(covars, collapse=" + "))), data = df[keep, , drop=FALSE])
  fit <- lm.fit(X, value[keep])
  out[keep] <- fit$residuals
  out
}

make_xy <- function(df, x_raw, y_raw, model, mode=c("overall","within_group")) {
  mode <- match.arg(mode)
  cov_primary <- c("age","bmi","disease_duration","constipation")
  cov_sens <- c(cov_primary, "levodopa_dose")

  if (model == "unadjusted") {
    df$x_used <- x_raw
    df$y_used <- y_raw
    return(df)
  }
  covars <- if (model=="primary") cov_primary else cov_sens

  if (mode=="overall") {
    df$x_used <- resid_on_covars(df, x_raw, covars)
    df$y_used <- resid_on_covars(df, y_raw, covars)
  } else {
    df$x_used <- NA_real_; df$y_used <- NA_real_
    for (g in levels(df$group)) {
      idx <- which(df$group==g)
      df$x_used[idx] <- resid_on_covars(df[idx, , drop=FALSE], x_raw[idx], covars)
      df$y_used[idx] <- resid_on_covars(df[idx, , drop=FALSE], y_raw[idx], covars)
    }
  }
  df
}

plot_groupwise_onefile <- function(df, title, xlab, ylab, out_pdf, out_png) {
  st_gr <- spearman_fast(df$x_used[df$group=="PD_GR"], df$y_used[df$group=="PD_GR"])
  st_mr <- spearman_fast(df$x_used[df$group=="PD_MR"], df$y_used[df$group=="PD_MR"])

  lab1 <- sprintf("GR: rho=%.3f, p=%.2g, n=%d", st_gr["rho"], st_gr["p"], as.integer(st_gr["n"]))
  lab2 <- sprintf("MR: rho=%.3f, p=%.2g, n=%d", st_mr["rho"], st_mr["p"], as.integer(st_mr["n"]))

  p <- ggplot(df, aes(x = x_used, y = y_used, color = group)) +
    geom_point(size = 2.6, alpha = 0.9) +
    geom_smooth(aes(group = group), method = "lm", se = TRUE, fill = "grey70", alpha = 0.25) +
    scale_color_manual(values = COL_MAP, breaks = names(LAB_MAP), labels = LAB_MAP) +
    theme_bw() +
    labs(title = title, x = xlab, y = ylab, color = "group") +
    annotate("text", x = Inf, y = Inf, label = lab1, hjust = 1.05, vjust = 2.0, size = 4) +
    annotate("text", x = Inf, y = Inf, label = lab2, hjust = 1.05, vjust = 3.2, size = 4)

  if (MAKE_PDF) ggsave(out_pdf, p, width = 7.8, height = 6.0)
  if (MAKE_PNG) ggsave(out_png, p, width = 7.8, height = 6.0, dpi = 300)
}

plot_overall_onefile <- function(df, title, xlab, ylab, out_pdf, out_png) {
  st <- spearman_fast(df$x_used, df$y_used)
  lab <- sprintf("Spearman rho=%.3f\np=%.2g\nn=%d", st["rho"], st["p"], as.integer(st["n"]))

  p <- ggplot(df, aes(x = x_used, y = y_used)) +
    geom_point(aes(color = group), size = 2.6, alpha = 0.9) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", fill = "grey70", alpha = 0.25) +
    scale_color_manual(values = COL_MAP, breaks = names(LAB_MAP), labels = LAB_MAP) +
    theme_bw() +
    labs(title = title, x = xlab, y = ylab, color = "group") +
    annotate("text", x = Inf, y = Inf, label = lab, hjust = 1.05, vjust = 1.2, size = 4)

  if (MAKE_PDF) ggsave(out_pdf, p, width = 7.8, height = 6.0)
  if (MAKE_PNG) ggsave(out_png, p, width = 7.8, height = 6.0, dpi = 300)
}

run_scatter_all_metrics <- function(df_alpha, metrics, dataset_tag, out_root) {
  dir_init(out_root)
  write_tsv0(meta, file.path(out_root, "00_manifest", "metadata_model_40.tsv"))
  write_tsv0(tibble(metric=metrics), file.path(out_root, "00_manifest", "alpha_metrics.tsv"))

  # keep needed columns
  df <- df_alpha |>
    select(sample_id, group, age, bmi, disease_duration, constipation, levodopa_dose, updrs3_reduction_pct, all_of(metrics)) |>
    filter(sample_id %in% sample_order) |>
    arrange(group, sample_id)

  for (mname in MODELS) {
    out_dir <- file.path(out_root, mname)
    dir_init(out_dir)

    for (metric in metrics) {
      x_raw <- df[[metric]]
      y_raw <- df$updrs3_reduction_pct

      # overall
      d_over <- make_xy(df, x_raw, y_raw, model=mname, mode="overall")
      write_tsv0(
        d_over |> select(sample_id, group, x_raw = all_of(metric), y_raw = updrs3_reduction_pct, x_used, y_used),
        file.path(out_dir, "per_sample", paste0("scatter_overall_data_", metric, ".tsv"))
      )
      plot_overall_onefile(
        d_over,
        title = paste0(dataset_tag, " ", metric, " vs UPDRS-III reduction (overall, ", mname, ")"),
        xlab = paste0(metric, " (adj residual, ", mname, ")"),
        ylab = paste0("UPDRS-III reduction (adj residual, ", mname, ")"),
        out_pdf = file.path(out_dir, "figures_pdf", paste0("scatter_overall_", metric, ".pdf")),
        out_png = file.path(out_dir, "figures_png", paste0("scatter_overall_", metric, ".png"))
      )

      # groupwise
      d_gw <- make_xy(df, x_raw, y_raw, model=mname, mode="within_group")
      write_tsv0(
        d_gw |> select(sample_id, group, x_raw = all_of(metric), y_raw = updrs3_reduction_pct, x_used, y_used),
        file.path(out_dir, "per_sample", paste0("scatter_groupwise_data_", metric, ".tsv"))
      )
      plot_groupwise_onefile(
        d_gw,
        title = paste0(dataset_tag, " ", metric, " vs UPDRS-III reduction (", mname, ")"),
        xlab = paste0(metric, " (adj residual, ", mname, ")"),
        ylab = paste0("UPDRS-III reduction (adj residual, ", mname, ")"),
        out_pdf = file.path(out_dir, "figures_pdf", paste0("scatter_groupwise_", metric, ".pdf")),
        out_png = file.path(out_dir, "figures_png", paste0("scatter_groupwise_", metric, ".png"))
      )
    }

    message("Done scatters: ", dataset_tag, " | ", mname)
  }
}

## -------------------------------
## 7) RUN
## -------------------------------
dir_init(OUT_ROOT)
write_tsv0(tibble(
  INPUT_16S_ALPHA = INPUT_16S_ALPHA,
  INPUT_VIROME_ALPHA = INPUT_VIROME_ALPHA,
  META_FILE = META_FILE
), file.path(OUT_ROOT, "00_manifest", "inputs.tsv"))

# 7.1 correlation + heatmap
run_alpha_clinical(alpha16, METRICS_16S, "16S_Enterococcus_within_genus", file.path(OUT_ROOT, "A_corr_heatmap_16S"))
run_alpha_clinical(alphav, METRICS_VIR, "Virome_EnterococcusHostSpecies", file.path(OUT_ROOT, "B_corr_heatmap_virome"))

# 7.2 UPDRS scatter (all metrics)
run_scatter_all_metrics(alpha16, METRICS_16S, "16S_Enterococcus_within_genus", file.path(OUT_ROOT, "C_scatter_UPDRS_16S"))
run_scatter_all_metrics(alphav, METRICS_VIR, "Virome_EnterococcusHostSpecies", file.path(OUT_ROOT, "D_scatter_UPDRS_virome"))

# session info
sink(file.path(OUT_ROOT, "00_manifest", "sessionInfo.txt"))
print(sessionInfo())
sink()

message("All done -> ", OUT_ROOT)