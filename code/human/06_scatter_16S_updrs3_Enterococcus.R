## =========================
## 09_scatter_16S_updrs3_Enterococcus.R
## One-figure-per-file scatter plots (PDF required)
## - X: CLR (or adjusted residual) of feature
## - Y: updrs3_reduction_pct (or adjusted residual)
## - Plot types per feature per model:
##   1) groupwise: 2 lines (PD_GR & PD_MR) on SAME plot
##   2) overall: 1 line (all samples) on SAME plot
## - Levels:
##   Genus: Enterococcus
##   Species: all species under genus Enterococcus
## - Models:
##   unadjusted / primary / sensitivity_add_dose
## Inputs in wd:
##   metadata.csv
##   table_filtered_w_tax.xls  或  table_f_even_depth_w_tax.xls
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

## ---------- params ----------
META_FILE <- "metadata.csv"
INPUT_16S <- "table_f_even_depth_w_tax.xls"  #抽平/未抽平

OUT_ROOT <- "results_16S_scatter_updrs3_even" #抽平/未抽平
dir.create(OUT_ROOT, showWarnings = FALSE)

GROUP_LEVELS <- c("PD_GR","PD_MR")  # GR first
PSEUDOCOUNT <- 0.5

MAKE_PDF <- TRUE
MAKE_PNG <- TRUE

# Colors consistent: GR green, MR red
COL_MAP <- c(PD_GR = "#7CAE00", PD_MR = "#F8766D")
LAB_MAP <- c(PD_GR = "GR", PD_MR = "MR")

MODELS <- c("unadjusted","primary","sensitivity_add_dose")

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
  dir.create(file.path(path, "figures_pdf"), showWarnings = FALSE)
  dir.create(file.path(path, "figures_png"), showWarnings = FALSE)
  dir.create(file.path(path, "logs"), showWarnings = FALSE)
}
write_tsv0 <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(x, path)
}
fill_unknown <- function(x, label) {
  x <- as.character(x)
  x[is.na(x) | x=="" | x=="NA"] <- paste0("Unassigned_", label)
  x
}
parse_percent <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "%", "")
  suppressWarnings(as.numeric(x))
}
safe_slug <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "[^A-Za-z0-9]+", "_")
  x <- str_replace_all(x, "_+", "_")
  x <- str_replace_all(x, "^_|_$", "")
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
    genus   = apply(parts, 1, get_rank, key = "g"),
    species = apply(parts, 1, get_rank, key = "s")
  )
}
aggregate_counts <- function(mat_feat_by_samp, feature_to_group_named) {
  grp <- feature_to_group_named[rownames(mat_feat_by_samp)]
  rowsum(mat_feat_by_samp, group = grp, reorder = TRUE)
}
clr_transform <- function(mat_feat_by_samp, pseudocount = 0.5) {
  x <- mat_feat_by_samp + pseudocount
  x <- sweep(x, 2, colSums(x), "/")
  lx <- log(x)
  sweep(lx, 2, colMeans(lx), "-")
}

# residualize vector on covariates (no group term)
resid_on_covars <- function(df, value, covars) {
  keep <- complete.cases(df[, covars, drop = FALSE]) & is.finite(value)
  out <- rep(NA_real_, length(value))
  if (sum(keep) < 6) return(out)
  X <- model.matrix(as.formula(paste("~", paste(covars, collapse = " + "))), data = df[keep, , drop = FALSE])
  fit <- lm.fit(X, value[keep])
  out[keep] <- fit$residuals
  out
}

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

make_xy <- function(df, x_raw, y_raw, model, mode = c("overall","within_group")) {
  mode <- match.arg(mode)

  cov_primary <- c("age","bmi","disease_duration","constipation")
  cov_sens <- c(cov_primary, "levodopa_dose")

  if (model == "unadjusted") {
    df$x_used <- x_raw
    df$y_used <- y_raw
    return(df)
  }

  covars <- if (model == "primary") cov_primary else cov_sens

  if (mode == "overall") {
    df$x_used <- resid_on_covars(df, x_raw, covars)
    df$y_used <- resid_on_covars(df, y_raw, covars)
  } else {
    df$x_used <- NA_real_
    df$y_used <- NA_real_
    for (g in levels(df$group)) {
      idx <- which(df$group == g)
      df$x_used[idx] <- resid_on_covars(df[idx, , drop=FALSE], x_raw[idx], covars)
      df$y_used[idx] <- resid_on_covars(df[idx, , drop=FALSE], y_raw[idx], covars)
    }
  }
  df
}

plot_groupwise_onefile <- function(df, title, xlab, ylab, out_pdf, out_png) {
  # Spearman by group (on x_used/y_used)
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

## ---------- load metadata ----------
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
  mutate(sample_id = normalize_id(.data[[sample_col]]),
         group = .data[[group_col]]) |>
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

stopifnot(nrow(meta) == 40)
sample_order <- meta$sample_id

write_tsv0(meta, file.path(OUT_ROOT, "00_manifest_metadata_model_40.tsv"))

## ---------- read 16S table ----------
dt <- fread(INPUT_16S, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
id_col  <- if ("#ID" %in% names(dt)) "#ID" else names(dt)[1]
tax_col <- if ("taxonomy" %in% names(dt)) "taxonomy" else if ("taxon" %in% names(dt)) "taxon" else NA_character_
stopifnot(!is.na(tax_col))

sample_cols <- names(dt)[str_detect(names(dt), "^PD\\d+$")]
keep_samples <- intersect(sample_order, sample_cols)
stopifnot(length(keep_samples) == 40)
keep_samples <- sample_order

asv_id <- as.character(dt[[id_col]])
mat_asv <- dt[, keep_samples, drop=FALSE]
mat_asv[] <- lapply(mat_asv, as.numeric)
mat_asv <- as.matrix(mat_asv)
rownames(mat_asv) <- asv_id

tax_tbl <- tibble(asv_id = asv_id, taxonomy_raw = dt[[tax_col]]) |>
  bind_cols(parse_tax_16s(dt[[tax_col]])) |>
  mutate(
    genus = fill_unknown(genus, "genus"),
    species = fill_unknown(species, "species")
  )

## ---------- build genus/species matrices + CLR ----------
mat_genus <- aggregate_counts(mat_asv, setNames(tax_tbl$genus, tax_tbl$asv_id))
mat_species_all <- aggregate_counts(mat_asv, setNames(tax_tbl$species, tax_tbl$asv_id))

clr_genus <- clr_transform(mat_genus, pseudocount = PSEUDOCOUNT)
clr_species <- clr_transform(mat_species_all, pseudocount = PSEUDOCOUNT)

if (!("Enterococcus" %in% rownames(clr_genus))) stop("Genus 'Enterococcus' not found in genus matrix.")

enterococcus_species <- rownames(clr_species)[str_detect(rownames(clr_species), "^Enterococcus\\b")]
if (length(enterococcus_species) == 0) stop("No Enterococcus species found (species names not starting with 'Enterococcus').")

## ---------- output folders ----------
genus_root <- file.path(OUT_ROOT, "Genus_Enterococcus")
species_root <- file.path(OUT_ROOT, "Species_Enterococcus")
dir_init(genus_root); dir_init(species_root)

write_tsv0(tibble(species = enterococcus_species),
           file.path(species_root, "tables", "enterococcus_species_list.tsv"))

## ---------- run all models ----------
for (m in MODELS) {

  ## ===== Genus Enterococcus =====
  out_dir <- file.path(genus_root, m)
  dir_init(out_dir)

  x_raw <- as.numeric(clr_genus["Enterococcus", sample_order])
  y_raw <- meta$updrs3_reduction_pct

  df_g <- make_xy(meta, x_raw, y_raw, model = m, mode = "within_group")
  df_o <- make_xy(meta, x_raw, y_raw, model = m, mode = "overall")

  write_tsv0(df_g |> mutate(feature="Enterococcus", level="Genus", model=m, plot_type="groupwise"),
             file.path(out_dir, "tables", "scatter_data_groupwise.tsv"))
  write_tsv0(df_o |> mutate(feature="Enterococcus", level="Genus", model=m, plot_type="overall"),
             file.path(out_dir, "tables", "scatter_data_overall.tsv"))

  xlab <- if (m=="unadjusted") "Enterococcus (CLR, genus level)" else paste0("Enterococcus (adj residual, ", m, ")")
  ylab <- if (m=="unadjusted") "MDS-UPDRS-III reduction (%)" else paste0("UPDRS-III reduction (adj residual, ", m, ")")

  plot_groupwise_onefile(df_g,
    title = paste0("Genus Enterococcus vs UPDRS-III reduction (", m, ")"),
    xlab=xlab, ylab=ylab,
    out_pdf=file.path(out_dir, "figures_pdf", "Genus_Enterococcus_groupwise.pdf"),
    out_png=file.path(out_dir, "figures_png", "Genus_Enterococcus_groupwise.png"))

  plot_overall_onefile(df_o,
    title = paste0("Genus Enterococcus vs UPDRS-III reduction (overall, ", m, ")"),
    xlab=xlab, ylab=ylab,
    out_pdf=file.path(out_dir, "figures_pdf", "Genus_Enterococcus_overall.pdf"),
    out_png=file.path(out_dir, "figures_png", "Genus_Enterococcus_overall.png"))

  ## ===== Species: all Enterococcus species =====
  out_dir2 <- file.path(species_root, m)
  dir_init(out_dir2)

  sum_rows <- list()
  long_g_all <- list()
  long_o_all <- list()

  for (sp in enterococcus_species) {
    sp_slug <- safe_slug(sp)

    x_raw_sp <- as.numeric(clr_species[sp, sample_order])
    y_raw_sp <- meta$updrs3_reduction_pct

    df_gs <- make_xy(meta, x_raw_sp, y_raw_sp, model = m, mode = "within_group") |>
      mutate(feature = sp, level="Species", model=m, plot_type="groupwise")
    df_os <- make_xy(meta, x_raw_sp, y_raw_sp, model = m, mode = "overall") |>
      mutate(feature = sp, level="Species", model=m, plot_type="overall")

    long_g_all[[length(long_g_all)+1]] <- df_gs
    long_o_all[[length(long_o_all)+1]] <- df_os

    st_gr <- spearman_fast(df_gs$x_used[df_gs$group=="PD_GR"], df_gs$y_used[df_gs$group=="PD_GR"])
    st_mr <- spearman_fast(df_gs$x_used[df_gs$group=="PD_MR"], df_gs$y_used[df_gs$group=="PD_MR"])
    st_all <- spearman_fast(df_os$x_used, df_os$y_used)

    sum_rows[[length(sum_rows)+1]] <- tibble(
      feature = sp, model = m,
      rho_GR = st_gr["rho"], p_GR = st_gr["p"], n_GR = st_gr["n"],
      rho_MR = st_mr["rho"], p_MR = st_mr["p"], n_MR = st_mr["n"],
      rho_all = st_all["rho"], p_all = st_all["p"], n_all = st_all["n"]
    )

    xlab_sp <- if (m=="unadjusted") paste0(sp, " (CLR, species)") else paste0(sp, " (adj residual, ", m, ")")
    ylab_sp <- if (m=="unadjusted") "MDS-UPDRS-III reduction (%)" else paste0("UPDRS-III reduction (adj residual, ", m, ")")

    # One file per plot
    plot_groupwise_onefile(df_gs,
      title = paste0("Species ", sp, " vs UPDRS-III reduction (", m, ")"),
      xlab=xlab_sp, ylab=ylab_sp,
      out_pdf=file.path(out_dir2, "figures_pdf", paste0("Species_", sp_slug, "_groupwise.pdf")),
      out_png=file.path(out_dir2, "figures_png", paste0("Species_", sp_slug, "_groupwise.png")))

    plot_overall_onefile(df_os,
      title = paste0("Species ", sp, " vs UPDRS-III reduction (overall, ", m, ")"),
      xlab=xlab_sp, ylab=ylab_sp,
      out_pdf=file.path(out_dir2, "figures_pdf", paste0("Species_", sp_slug, "_overall.pdf")),
      out_png=file.path(out_dir2, "figures_png", paste0("Species_", sp_slug, "_overall.png")))
  }

  # Save per-sample long tables and summary (per model)
  write_tsv0(bind_rows(long_g_all),
             file.path(out_dir2, "tables", "scatter_data_groupwise_all_species_long.tsv"))
  write_tsv0(bind_rows(long_o_all),
             file.path(out_dir2, "tables", "scatter_data_overall_all_species_long.tsv"))
  write_tsv0(bind_rows(sum_rows),
             file.path(out_dir2, "tables", "species_correlation_summary.tsv"))

  message("Done model: ", m)
}

message("All scatterplots done -> ", OUT_ROOT)