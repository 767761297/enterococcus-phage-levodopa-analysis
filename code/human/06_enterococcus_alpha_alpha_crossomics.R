#!/usr/bin/env Rscript
## =========================================================
## Cross-omics alpha-alpha correlation
## X = Enterococcus-host phage alpha (virome; Host_species level)
## Y = Enterococcus alpha (16S; within-genus)
##
## Outputs per model:
## - tables: rho/p/q matrices + long table (all pairs)
## - figures: heatmap (X=virome metrics, Y=16S metrics) with q-significance marks
## - figures: scatter (overall + groupwise) for ALL metrics or MATCHED pairs
##
## Models:
## - unadjusted: raw metrics
## - primary: residualize BOTH sides on group + age + bmi + disease_duration + constipation
## - sensitivity_add_dose: primary + levodopa_dose
##
## Notes:
## - For groupwise scatter, residualization is done WITHIN group and WITHOUT group term (same spirit as your old scatter script).
## - Robust to NA/NaN/Inf: residualization is per-metric with complete-case per metric.
## =========================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(janitor)
  library(ggplot2)
  library(pheatmap)
})

## ------------------------
## 0) USER SETTINGS
## ------------------------
# Per-sample alpha tables (generated earlier)
INPUT_16S_ALPHA    <- "alpha_per_sample_Enterococcus_species.tsv"
INPUT_VIROME_ALPHA <- "alpha_per_sample_EnterococcusHostSpeciesPhage.tsv"
META_FILE          <- "metadata.csv"

# Output root folder
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_ROOT <- paste0("results_alpha_alpha_crossomics_", timestamp)

# Groups and colors (match your convention)
GROUP_LEVELS <- c("PD_GR", "PD_MR")
COL_MAP <- c(PD_GR = "#7CAE00", PD_MR = "#F8766D")  # GR green, MR red
LAB_MAP <- c(PD_GR = "GR", PD_MR = "MR")

MODELS <- c("unadjusted", "primary", "sensitivity_add_dose")

# Heatmap significance markers based on q (global BH across all pairs)
Q_CUTS <- c(0.05, 0.10)   # ** for <0.05; * for <0.10

# Scatter plotting:
# "matched" = only matched metric pairs (recommended)
# "all"     = all combinations (may generate many plots)
SCATTER_MODE <- "matched"

MAKE_PDF <- TRUE
MAKE_PNG <- TRUE

## ------------------------
## 1) Helpers
## ------------------------
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

# residualize numeric vector on design matrix (lm.fit), per-variable NA-safe
residualize_vector <- function(y, X) {
  ok <- complete.cases(X) & is.finite(y)
  out <- rep(NA_real_, length(y))
  if (sum(ok) < 6) return(out)
  fit <- lm.fit(X[ok, , drop=FALSE], y[ok])
  out[ok] <- fit$residuals
  out
}

# Spearman rho + t approximation p (NA-safe)
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

# Build adjustment design matrix for OVERALL residualization (includes group)
build_X_overall <- function(meta, model_name, remove_terms = character()) {
  if (model_name == "unadjusted") return(NULL)
  base <- c("group", "age", "bmi", "disease_duration", "constipation")
  if (model_name == "sensitivity_add_dose") base <- c(base, "levodopa_dose")
  base <- setdiff(base, remove_terms)
  f <- as.formula(paste("~", paste(base, collapse=" + ")))
  model.matrix(f, data = meta)
}

# Build adjustment design matrix for WITHIN-GROUP residualization (no group term)
build_X_within <- function(meta_g, model_name, remove_terms = character()) {
  if (model_name == "unadjusted") return(NULL)
  base <- c("age", "bmi", "disease_duration", "constipation")
  if (model_name == "sensitivity_add_dose") base <- c(base, "levodopa_dose")
  base <- setdiff(base, remove_terms)
  f <- as.formula(paste("~", paste(base, collapse=" + ")))
  model.matrix(f, data = meta_g)
}

## ------------------------
## 2) Load metadata (robust column mapping)
## ------------------------
meta_raw <- read_csv_smart(META_FILE) |> janitor::clean_names()

sample_col <- find_col(meta_raw, c("^sample$"))
group_col  <- find_col(meta_raw, c("^group$"))
age_col    <- find_col(meta_raw, c("^age$"))
bmi_col    <- find_col(meta_raw, c("^bmi"))
dur_col    <- find_col(meta_raw, c("disease_duration"))
con_col    <- find_col(meta_raw, c("constipation"))
dose_col   <- find_col(meta_raw, c("levodopa.*dose|levodopa_daily_dose"))
updrs_col  <- find_col(meta_raw, c("mds.*updrs.*iii.*reduction", "updrs.*iii.*reduction", "reduction.*updrs.*iii"))

stopifnot(!is.na(sample_col), !is.na(group_col), !is.na(age_col),
          !is.na(bmi_col), !is.na(dur_col), !is.na(con_col), !is.na(dose_col))

meta <- meta_raw |>
  mutate(
    sample_id = normalize_id(.data[[sample_col]]),
    group = factor(.data[[group_col]], levels = GROUP_LEVELS),
    age = suppressWarnings(as.numeric(.data[[age_col]])),
    bmi = suppressWarnings(as.numeric(.data[[bmi_col]])),
    disease_duration = suppressWarnings(as.numeric(.data[[dur_col]])),
    constipation = factor(as.integer(.data[[con_col]]), levels=c(0,1), labels=c("No","Yes")),
    levodopa_dose = suppressWarnings(as.numeric(.data[[dose_col]])),
    updrs3_reduction_pct = if (!is.na(updrs_col)) parse_percent(.data[[updrs_col]]) else NA_real_
  ) |>
  filter(group %in% GROUP_LEVELS) |>
  arrange(group, sample_id)

sample_order <- meta$sample_id

## ------------------------
## 3) Load alpha tables (merge metadata if needed)
## ------------------------
load_alpha <- function(path, dataset_name) {
  df <- read_tsv_smart(path) |> janitor::clean_names()
  sid <- find_col(df, c("^sample_id$", "^sample$"))
  stopifnot(!is.na(sid))
  df <- df |> mutate(sample_id = normalize_id(.data[[sid]]))

  need_merge <- any(!c("group","age","bmi","disease_duration","constipation","levodopa_dose") %in% names(df))
  if (need_merge) {
    df <- df |> select(-matches("^group$|^age$|^bmi$|^disease_duration$|^constipation$|^levodopa_dose$|^updrs3_reduction"))
    df <- df |> left_join(meta, by="sample_id")
  } else {
    df <- df |>
      mutate(
        group = factor(as.character(group), levels = GROUP_LEVELS),
        age = suppressWarnings(as.numeric(age)),
        bmi = suppressWarnings(as.numeric(bmi)),
        disease_duration = suppressWarnings(as.numeric(disease_duration)),
        constipation = factor(as.character(constipation), levels=c("No","Yes")),
        levodopa_dose = suppressWarnings(as.numeric(levodopa_dose))
      )
  }

  df <- df |> filter(sample_id %in% sample_order) |> arrange(group, sample_id)
  df$dataset <- dataset_name
  df
}

alpha16 <- load_alpha(INPUT_16S_ALPHA, "16S")
alphav  <- load_alpha(INPUT_VIROME_ALPHA, "Virome")

## Metrics (you can edit if needed)
METRICS_16S <- c("observed","chao1","shannon","exp_shannon","simpson","inv_simpson","pielou_evenness")
METRICS_VIR <- c("observed_gt0","observed_gt1","shannon","exp_shannon","simpson","inv_simpson","pielou_evenness")

# sanity: keep only existing columns
METRICS_16S <- METRICS_16S[METRICS_16S %in% names(alpha16)]
METRICS_VIR <- METRICS_VIR[METRICS_VIR %in% names(alphav)]

stopifnot(length(METRICS_16S) > 0, length(METRICS_VIR) > 0)

## Optional: prefix for clarity (used in heatmap/scatter labels)
METRICS_16S_LAB <- paste0("16S_", METRICS_16S)
METRICS_VIR_LAB <- paste0("VIR_", METRICS_VIR)

## Build matched pairs (for scatter)
matched_pairs <- tibble(
  metric_16S = c("observed","chao1","shannon","exp_shannon","simpson","inv_simpson","pielou_evenness"),
  metric_vir = c("observed_gt1","observed_gt1","shannon","exp_shannon","simpson","inv_simpson","pielou_evenness")
) |>
  filter(metric_16S %in% METRICS_16S, metric_vir %in% METRICS_VIR)

## ------------------------
## 4) Core: alpha-alpha correlation per model
## ------------------------
run_alpha_alpha <- function(model_name, out_dir) {
  dir_init(out_dir)

  # matrices aligned to sample_order
  M16 <- alpha16 |>
    select(sample_id, all_of(METRICS_16S)) |>
    arrange(match(sample_id, sample_order)) |>
    column_to_rownames("sample_id") |>
    as.matrix()
  MV <- alphav |>
    select(sample_id, all_of(METRICS_VIR)) |>
    arrange(match(sample_id, sample_order)) |>
    column_to_rownames("sample_id") |>
    as.matrix()

  # Clean non-finite to NA
  M16[!is.finite(M16)] <- NA_real_
  MV[!is.finite(MV)] <- NA_real_

  # Residualization (OVERALL; includes group, like your old clinical corr script)
  if (model_name != "unadjusted") {
    X <- build_X_overall(meta, model_name)
    # residualize each metric separately (NA-safe)
    for (j in seq_len(ncol(M16))) M16[, j] <- residualize_vector(M16[, j], X)
    for (j in seq_len(ncol(MV)))  MV[, j]  <- residualize_vector(MV[, j],  X)
  }

  # Compute pairwise correlations: Y=16S metrics (rows), X=VIR metrics (cols)
  rho <- matrix(NA_real_, nrow=length(METRICS_16S), ncol=length(METRICS_VIR),
                dimnames=list(METRICS_16S_LAB, METRICS_VIR_LAB))
  pmat <- rho
  nmat <- rho

  long_list <- list()
  idx <- 0L

  for (i in seq_along(METRICS_16S)) {
    for (j in seq_along(METRICS_VIR)) {
      x <- MV[, j]
      y <- M16[, i]
      st <- spearman_fast(x, y)
      idx <- idx + 1L
      long_list[[idx]] <- tibble(
        model = model_name,
        metric_16S = METRICS_16S[i],
        metric_virome = METRICS_VIR[j],
        metric_16S_lab = METRICS_16S_LAB[i],
        metric_virome_lab = METRICS_VIR_LAB[j],
        rho = as.numeric(st["rho"]),
        p_value = as.numeric(st["p"]),
        n_used = as.integer(st["n"])
      )
      rho[METRICS_16S_LAB[i], METRICS_VIR_LAB[j]] <- as.numeric(st["rho"])
      pmat[METRICS_16S_LAB[i], METRICS_VIR_LAB[j]] <- as.numeric(st["p"])
      nmat[METRICS_16S_LAB[i], METRICS_VIR_LAB[j]] <- as.numeric(st["n"])
    }
  }

  long_df <- bind_rows(long_list)
  # Global BH across ALL pairs (recommended for this matrix)
  long_df <- long_df |> mutate(q_value = p.adjust(p_value, method="BH"))

  # q matrix
  qmat <- matrix(NA_real_, nrow=nrow(rho), ncol=ncol(rho), dimnames=dimnames(rho))
  for (k in seq_len(nrow(long_df))) {
    qmat[ long_df$metric_16S_lab[k], long_df$metric_virome_lab[k] ] <- long_df$q_value[k]
  }

  # Save tables
  write_tsv0(long_df, file.path(out_dir, "tables", "alpha_alpha_long.tsv"))
  write_tsv0(as.data.frame(rho)  |> rownames_to_column("metric_16S"),
             file.path(out_dir, "tables", "rho_matrix.tsv"))
  write_tsv0(as.data.frame(pmat) |> rownames_to_column("metric_16S"),
             file.path(out_dir, "tables", "p_matrix.tsv"))
  write_tsv0(as.data.frame(qmat) |> rownames_to_column("metric_16S"),
             file.path(out_dir, "tables", "q_matrix.tsv"))
  write_tsv0(as.data.frame(nmat) |> rownames_to_column("metric_16S"),
             file.path(out_dir, "tables", "n_matrix.tsv"))

  # Heatmap: display_numbers by q cutoffs
  sig <- matrix("", nrow=nrow(qmat), ncol=ncol(qmat), dimnames=dimnames(qmat))
  sig[qmat < Q_CUTS[2]] <- "*"
  sig[qmat < Q_CUTS[1]] <- "**"

  hm_title <- paste0("Alpha–alpha coupling | ", model_name,
                     " | X=Enterococcus-host phage (virome), Y=Enterococcus (16S)")

  if (MAKE_PDF) {
    pdf(file.path(out_dir, "figures_pdf", "heatmap_rho.pdf"),
        width = 8.5, height = max(4.0, 0.35*nrow(rho) + 2.5))
    pheatmap(
      rho,
      cluster_rows = FALSE, cluster_cols = FALSE,
      display_numbers = sig, number_color = "black",
      main = hm_title, fontsize = 11
    )
    dev.off()
  }
  if (MAKE_PNG) {
    png(file.path(out_dir, "figures_png", "heatmap_rho.png"),
        width = 1400, height = max(800, 150*nrow(rho)), res = 150)
    pheatmap(
      rho,
      cluster_rows = FALSE, cluster_cols = FALSE,
      display_numbers = sig, number_color = "black",
      main = hm_title, fontsize = 11
    )
    dev.off()
  }

  invisible(list(long=long_df, rho=rho, p=pmat, q=qmat))
}

## ------------------------
## 5) Scatter plots (overall + groupwise) per model
## ------------------------
# residualize within-group (no group term) for scatter (to mimic your old scatter behavior)
resid_scatter <- function(df, x_raw, y_raw, model_name, mode=c("overall","within_group")) {
  mode <- match.arg(mode)

  if (model_name == "unadjusted") {
    df$x_used <- x_raw
    df$y_used <- y_raw
    return(df)
  }

  cov_primary <- c("age","bmi","disease_duration","constipation")
  cov_sens <- c(cov_primary, "levodopa_dose")
  covars <- if (model_name == "primary") cov_primary else cov_sens

  if (mode == "overall") {
    X <- model.matrix(as.formula(paste("~", paste(covars, collapse=" + "))), data=df)
    df$x_used <- residualize_vector(x_raw, X)
    df$y_used <- residualize_vector(y_raw, X)
  } else {
    df$x_used <- NA_real_; df$y_used <- NA_real_
    for (g in levels(df$group)) {
      idx <- which(df$group == g)
      if (length(idx) < 6) next
      dfg <- df[idx, , drop=FALSE]
      Xg <- model.matrix(as.formula(paste("~", paste(covars, collapse=" + "))), data=dfg)
      df$x_used[idx] <- residualize_vector(x_raw[idx], Xg)
      df$y_used[idx] <- residualize_vector(y_raw[idx], Xg)
    }
  }
  df
}

plot_overall <- function(df, title, xlab, ylab, out_pdf, out_png) {
  df <- df |> filter(is.finite(x_used), is.finite(y_used))
  st <- spearman_fast(df$x_used, df$y_used)
  lab <- sprintf("Spearman rho=%.3f\np=%.2g\nn=%d",
                 st["rho"], st["p"], as.integer(st["n"]))

  p <- ggplot(df, aes(x=x_used, y=y_used)) +
    geom_point(aes(color=group), size=2.6, alpha=0.9) +
    geom_smooth(aes(group=1), method="lm", se=TRUE,
                color="black", fill="grey70", alpha=0.25) +
    scale_color_manual(values=COL_MAP, breaks=names(LAB_MAP), labels=LAB_MAP) +
    theme_bw() +
    labs(title=title, x=xlab, y=ylab, color="group") +
    annotate("text", x=Inf, y=Inf, label=lab, hjust=1.05, vjust=1.2, size=4)

  if (MAKE_PDF) ggsave(out_pdf, p, width=7.8, height=6.0)
  if (MAKE_PNG) ggsave(out_png, p, width=7.8, height=6.0, dpi=300)
}

plot_groupwise <- function(df, title, xlab, ylab, out_pdf, out_png) {
  df <- df |> filter(is.finite(x_used), is.finite(y_used))
  st_gr <- spearman_fast(df$x_used[df$group=="PD_GR"], df$y_used[df$group=="PD_GR"])
  st_mr <- spearman_fast(df$x_used[df$group=="PD_MR"], df$y_used[df$group=="PD_MR"])
  lab1 <- sprintf("GR: rho=%.3f, p=%.2g, n=%d", st_gr["rho"], st_gr["p"], as.integer(st_gr["n"]))
  lab2 <- sprintf("MR: rho=%.3f, p=%.2g, n=%d", st_mr["rho"], st_mr["p"], as.integer(st_mr["n"]))

  p <- ggplot(df, aes(x=x_used, y=y_used, color=group)) +
    geom_point(size=2.6, alpha=0.9) +
    geom_smooth(aes(group=group), method="lm", se=TRUE, fill="grey70", alpha=0.25) +
    scale_color_manual(values=COL_MAP, breaks=names(LAB_MAP), labels=LAB_MAP) +
    theme_bw() +
    labs(title=title, x=xlab, y=ylab, color="group") +
    annotate("text", x=Inf, y=Inf, label=lab1, hjust=1.05, vjust=2.0, size=4) +
    annotate("text", x=Inf, y=Inf, label=lab2, hjust=1.05, vjust=3.2, size=4)

  if (MAKE_PDF) ggsave(out_pdf, p, width=7.8, height=6.0)
  if (MAKE_PNG) ggsave(out_png, p, width=7.8, height=6.0, dpi=300)
}

run_scatters <- function(model_name, out_dir) {
  dir_init(out_dir)

  df <- meta |>
    select(sample_id, group, age, bmi, disease_duration, constipation, levodopa_dose) |>
    left_join(
      alpha16 |>
        select(sample_id, all_of(METRICS_16S)) |>
        rename_with(~ paste0("16S_", .x), all_of(METRICS_16S)),
      by = "sample_id"
    ) |>
    left_join(
      alphav |>
        select(sample_id, all_of(METRICS_VIR)) |>
        rename_with(~ paste0("VIR_", .x), all_of(METRICS_VIR)),
      by = "sample_id"
    ) |>
    arrange(group, sample_id)

  # choose pairs
  pairs <- if (SCATTER_MODE == "all") {
    expand_grid(
      metric_16S = paste0("16S_", METRICS_16S),
      metric_vir = paste0("VIR_", METRICS_VIR)
    )
  } else {
    matched_pairs |>
      transmute(
        metric_16S = paste0("16S_", metric_16S),
        metric_vir = paste0("VIR_", metric_vir)
      )
  }

  # loop
  for (k in seq_len(nrow(pairs))) {
    m16 <- pairs$metric_16S[k]
    mvr <- pairs$metric_vir[k]

    x_raw <- df[[mvr]]
    y_raw <- df[[m16]]

    # overall scatter
    d_over <- resid_scatter(df, x_raw, y_raw, model_name, mode="overall")
    write_tsv0(
      d_over |> transmute(sample_id, group,
                          x_raw = .data[[mvr]], y_raw = .data[[m16]],
                          x_used, y_used),
      file.path(out_dir, "per_sample", paste0("scatter_overall_", "VIR_", mvr, "__", "16S_", m16, ".tsv"))
    )

    plot_overall(
      d_over,
      title = paste0("Overall | ", model_name, " | X=VIR_", mvr, " ; Y=16S_", m16),
      xlab  = paste0("VIR_", mvr, " (adj residual, ", model_name, ")"),
      ylab  = paste0("16S_", m16, " (adj residual, ", model_name, ")"),
      out_pdf = file.path(out_dir, "figures_pdf", paste0("scatter_overall__", "VIR_", mvr, "__", "16S_", m16, ".pdf")),
      out_png = file.path(out_dir, "figures_png", paste0("scatter_overall__", "VIR_", mvr, "__", "16S_", m16, ".png"))
    )

    # groupwise scatter
    d_gw <- resid_scatter(df, x_raw, y_raw, model_name, mode="within_group")
    write_tsv0(
      d_gw |> transmute(sample_id, group,
                        x_raw = .data[[mvr]], y_raw = .data[[m16]],
                        x_used, y_used),
      file.path(out_dir, "per_sample", paste0("scatter_groupwise_", "VIR_", mvr, "__", "16S_", m16, ".tsv"))
    )

    plot_groupwise(
      d_gw,
      title = paste0("Groupwise | ", model_name, " | X=VIR_", mvr, " ; Y=16S_", m16),
      xlab  = paste0("VIR_", mvr, " (adj residual, ", model_name, ")"),
      ylab  = paste0("16S_", m16, " (adj residual, ", model_name, ")"),
      out_pdf = file.path(out_dir, "figures_pdf", paste0("scatter_groupwise__", "VIR_", mvr, "__", "16S_", m16, ".pdf")),
      out_png = file.path(out_dir, "figures_png", paste0("scatter_groupwise__", "VIR_", mvr, "__", "16S_", m16, ".png"))
    )
  }
}

## ------------------------
## 6) RUN ALL
## ------------------------
dir_init(OUT_ROOT)

# manifest
write_tsv0(meta, file.path(OUT_ROOT, "00_manifest", "metadata_model_40.tsv"))
write_tsv0(tibble(
  INPUT_16S_ALPHA = INPUT_16S_ALPHA,
  INPUT_VIROME_ALPHA = INPUT_VIROME_ALPHA,
  META_FILE = META_FILE,
  SCATTER_MODE = SCATTER_MODE
), file.path(OUT_ROOT, "00_manifest", "inputs.tsv"))

write_tsv0(tibble(metric_16S = METRICS_16S, metric_16S_lab = METRICS_16S_LAB),
           file.path(OUT_ROOT, "00_manifest", "metrics_16S.tsv"))
write_tsv0(tibble(metric_virome = METRICS_VIR, metric_virome_lab = METRICS_VIR_LAB),
           file.path(OUT_ROOT, "00_manifest", "metrics_virome.tsv"))
write_tsv0(matched_pairs, file.path(OUT_ROOT, "00_manifest", "matched_pairs.tsv"))

# per model
for (m in MODELS) {
  message("Running model: ", m)

  out_m <- file.path(OUT_ROOT, paste0("A_alpha_alpha_", m))
  run_alpha_alpha(m, out_m)

  out_s <- file.path(OUT_ROOT, paste0("B_scatter_", m))
  run_scatters(m, out_s)
}

sink(file.path(OUT_ROOT, "00_manifest", "sessionInfo.txt"))
print(sessionInfo())
sink()

message("All done -> ", OUT_ROOT)