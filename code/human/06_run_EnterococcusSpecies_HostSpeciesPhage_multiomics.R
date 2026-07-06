#!/usr/bin/env Rscript
## =========================================================
## 06_run_EnterococcusSpecies_HostSpeciesPhage_multiomics_NOstack_FIXED.R
##
## What it does (NO stacked bars):
## 1) 16S: Extract Enterococcus species (from Genus/Species columns OR taxonomy string),
##    aggregate to species, filter (prev>=10%, total>=10), CLR (p=0.5).
## 2) Virome: vOTU TPM + CHERRY host prediction -> aggregate to Host_species,
##    build TWO Y sets:
##      - All Host_species (your “full virome” definition at Host_species level)
##      - Enterococcus-host subset (host_genus == Enterococcus)
##    filter (prev>=10%), CLR (p=0.01), name = "<host>_virus"
## 3) Cross-omics correlation (X=16S Enterococcus species, Y=Host_species_virus):
##    Models:
##      - unadjusted: CLR
##      - primary: residualize on group + age + bmi + disease_duration + constipation
##      - sensitivity_add_dose: add levodopa_dose
##    Output full rho/p/q matrices + long table + per-sample matrices used.
## 4) Heatmaps:
##    - All Host_species: select VIRUS columns by minP across bacteria rows, keep columns with minP<0.05;
##      if too many -> top 50 by minP; if none -> top 50 by minP.
##    - Enterococcus-host subset: NO threshold; plot ALL columns (after prevalence filter).
## 5) Scatterplots (same style as your old script):
##    - For each model: (overall + groupwise) for matched pairs:
##      intersection of 16S species names and virome host_species names (remove "_virus").
##
## Bug fixes vs old scripts:
## - Robust metadata column matching (clean_names + µ replacement)
## - Robust reading of ".xls" TSV via fread_smart (UTF-8/GB18030)
## - Residualization avoids qr.resid NA crash (complete-case by model; feature-wise lm.fit)
## - Avoid column name collisions by explicit prefixing in outputs
## - Always create dirs; always output PDF + PNG
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
})

MAKE_PDF <- TRUE
MAKE_PNG <- TRUE

## -------------------------
## 0) Inputs (auto-pick 16S table)
## -------------------------
META_FILE <- "metadata.csv"

INPUT_16S_CANDIDATES <- c("table_f_even_depth_w_tax.xls", "table_filtered_w_tax.xls")
INPUT_16S <- INPUT_16S_CANDIDATES[file.exists(INPUT_16S_CANDIDATES)][1]
if (is.na(INPUT_16S)) stop("No 16S table found: table_f_even_depth_w_tax.xls or table_filtered_w_tax.xls")

TPM_FILE <- "vOTU_abundance_matrix_TPM.tsv"
CHERRY_FILE <- "cherry_prediction.tsv"

## -------------------------
## 1) Global params
## -------------------------
GROUP_LEVELS <- c("PD_GR", "PD_MR")
MODELS <- c("unadjusted", "primary", "sensitivity_add_dose")

# CLR pseudocounts (match old script)
PSEUDOCOUNT_16S <- 0.5
PSEUDOCOUNT_VIROME <- 0.01

# Filtering
MIN_PREV_ALL <- 0.10         # >=4/40
MIN_TOTAL_16S <- 10
MIN_TOTAL_VIROME <- 0

# Heatmap column selection for All Host_species
P_CUTOFF <- 0.05
HEATMAP_MAX_COLS <- 50

# Output
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_ROOT <- paste0("EHP_", timestamp)

# Scatter
DO_SCATTER <- TRUE

# Colors (match your style)
COL_MAP <- c(PD_GR="#7CAE00", PD_MR="#F8766D")
LAB_MAP <- c(PD_GR="GR", PD_MR="MR")

## -------------------------
## 2) I/O helpers
## -------------------------
dir_init <- function(path) {
  dir.create(path, recursive=TRUE, showWarnings=FALSE)
  dir.create(file.path(path, "00_manifest"), showWarnings=FALSE)
  dir.create(file.path(path, "tables"), showWarnings=FALSE)
  dir.create(file.path(path, "per_sample"), showWarnings=FALSE)
  dir.create(file.path(path, "figures_pdf"), showWarnings=FALSE)
  dir.create(file.path(path, "figures_png"), showWarnings=FALSE)
  dir.create(file.path(path, "logs"), showWarnings=FALSE)
}

write_tsv0 <- function(x, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  readr::write_tsv(x, path)
}

fread_smart <- function(path, sep="\t") {
  # robust encoding for Chinese Windows
  for (enc in c("UTF-8", "GB18030")) {
    try({
      dt <- data.table::fread(path, sep=sep, encoding=enc, data.table=FALSE, fill=TRUE, quote="")
      return(dt)
    }, silent=TRUE)
  }
  # fallback
  dt <- read.delim(path, sep=sep, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, fileEncoding="GB18030")
  dt
}

read_csv_smart <- function(path) {
  tryCatch(
    readr::read_csv(path, show_col_types = FALSE, locale = locale(encoding="GB18030")),
    error=function(e) readr::read_csv(path, show_col_types = FALSE, locale = locale(encoding="UTF-8"))
  )
}

normalize_id <- function(x) x |> as.character() |> str_trim() |> str_replace_all("\\s+", "")

## -------------------------
## 3) Metadata standardization (robust)
## -------------------------
meta_raw <- read_csv_smart(META_FILE)

# clean names with µ replacement to silence warning
meta <- meta_raw |>
  janitor::clean_names(replace = c("\u03bc"="u", "\u00b5"="u"))

find_col <- function(nms, patterns) {
  for (p in patterns) {
    hit <- nms[str_detect(nms, regex(p, ignore_case=TRUE))]
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

nms <- names(meta)
col_sample <- find_col(nms, c("^sample$", "sample_id", "sampleid"))
col_group  <- find_col(nms, c("^group$"))
col_age    <- find_col(nms, c("^age$"))
col_bmi    <- find_col(nms, c("^bmi"))
col_ddur   <- find_col(nms, c("disease_duration"))
col_const  <- find_col(nms, c("constipation"))
col_dose   <- find_col(nms, c("levodopa.*dose", "levodopa_daily_dose"))

if (any(is.na(c(col_sample,col_group,col_age,col_bmi,col_ddur,col_const,col_dose)))) {
  stop("metadata.csv missing required columns after clean_names.\n",
       "Found columns: ", paste(nms, collapse=", "), "\n",
       "Need at least: Sample, Group, Age, BMI, Disease duration, Constipation, Levodopa daily dose")
}

meta_std <- meta |>
  transmute(
    sample_id = normalize_id(.data[[col_sample]]),
    group = as.character(.data[[col_group]]),
    age = suppressWarnings(as.numeric(.data[[col_age]])),
    bmi = suppressWarnings(as.numeric(.data[[col_bmi]])),
    disease_duration = suppressWarnings(as.numeric(.data[[col_ddur]])),
    constipation = as.integer(.data[[col_const]]),
    levodopa_dose = suppressWarnings(as.numeric(.data[[col_dose]]))
  ) |>
  filter(group %in% GROUP_LEVELS) |>
  mutate(
    group = factor(group, levels=GROUP_LEVELS),
    constipation = factor(constipation, levels=c(0,1), labels=c("No","Yes"))
  ) |>
  arrange(group, sample_id)

stopifnot(nrow(meta_std) > 0)

sample_order <- meta_std$sample_id

## -------------------------
## 4) Taxonomy parsing for 16S
## -------------------------
is_sample_col <- function(x) str_detect(x, "^PD\\d+$")

parse_tax_rank <- function(tax, rank_char) {
  if (is.na(tax)) return("")
  parts <- str_split(as.character(tax), ";", simplify=TRUE)
  parts <- str_trim(parts)
  pref1 <- paste0(rank_char, "_")
  pref2 <- paste0(rank_char, "__")
  hit <- parts[str_starts(parts, pref2)]
  if (length(hit) > 0) return(str_remove(hit[1], fixed(pref2)))
  hit <- parts[str_starts(parts, pref1)]
  if (length(hit) > 0) return(str_remove(hit[1], fixed(pref1)))
  ""
}

normalize_species_label <- function(genus, sp) {
  sp <- str_trim(as.character(sp))
  if (is.na(sp) || sp=="" || tolower(sp) %in% c("na","nan","unclassified","uncultured","unknown")) {
    return(paste0(genus, "_unclassified"))
  }
  sp <- str_remove(sp, "^s__")
  if (!str_detect(tolower(sp), "^enterococcus")) return(paste0(genus, " ", sp))
  sp
}

## CLR transform (features x samples)
clr_transform <- function(mat_feat_by_samp, pseudocount=0.5) {
  mat <- as.matrix(mat_feat_by_samp)
  storage.mode(mat) <- "numeric"
  mat <- mat + pseudocount
  cs <- colSums(mat)
  cs[cs <= 0] <- NA_real_
  ra <- sweep(mat, 2, cs, "/")
  ra[!is.finite(ra)] <- NA_real_
  logra <- log(ra)
  # center per sample (column)
  mu <- colMeans(logra, na.rm=TRUE)
  clr <- sweep(logra, 2, mu, "-")
  clr
}

## residualize samples x features matrix by design X (feature-wise, NA-safe)
residualize_matrix_safe <- function(Y_samp_by_feat, X_design) {
  stopifnot(nrow(Y_samp_by_feat) == nrow(X_design))
  out <- matrix(NA_real_, nrow=nrow(Y_samp_by_feat), ncol=ncol(Y_samp_by_feat),
                dimnames=dimnames(Y_samp_by_feat))
  okX <- complete.cases(X_design)
  if (sum(okX) < 6) return(out)
  for (j in seq_len(ncol(Y_samp_by_feat))) {
    y <- Y_samp_by_feat[, j]
    ok <- okX & is.finite(y)
    if (sum(ok) < 6) next
    fit <- lm.fit(X_design[ok, , drop=FALSE], y[ok])
    out[ok, j] <- fit$residuals
  }
  out
}

build_design <- function(meta_df, model_name) {
  if (model_name == "primary") {
    model.matrix(~ group + age + bmi + disease_duration + constipation, data=meta_df)
  } else if (model_name == "sensitivity_add_dose") {
    model.matrix(~ group + age + bmi + disease_duration + constipation + levodopa_dose, data=meta_df)
  } else {
    model.matrix(~ 1, data=meta_df)
  }
}

spearman_fast <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 6) return(c(rho=NA_real_, p=NA_real_, n=n))
  rho <- suppressWarnings(cor(x, y, method="spearman"))
  rho2 <- min(0.999999, max(-0.999999, rho))
  tval <- rho2 * sqrt((n - 2) / (1 - rho2^2))
  p <- 2 * pt(-abs(tval), df=n-2)
  c(rho=rho, p=p, n=n)
}

plot_heatmap <- function(mat, out_pdf, out_png, title=NULL) {
  dir.create(dirname(out_pdf), recursive=TRUE, showWarnings=FALSE)
  dir.create(dirname(out_png), recursive=TRUE, showWarnings=FALSE)
  if (requireNamespace("pheatmap", quietly=TRUE)) {
    cr <- (nrow(mat) >= 2)
    cc <- (ncol(mat) >= 2)
    if (MAKE_PDF) {
      pdf(out_pdf, width=10, height=max(4, 0.12*nrow(mat) + 2))
      pheatmap::pheatmap(mat, cluster_rows=cr, cluster_cols=cc, main=title, fontsize_row=7, fontsize_col=7)
      dev.off()
    }
    if (MAKE_PNG) {
      png(out_png, width=2600, height=max(1600, 28*nrow(mat) + 800), res=300)
      pheatmap::pheatmap(mat, cluster_rows=cr, cluster_cols=cc, main=title, fontsize_row=7, fontsize_col=7)
      dev.off()
    }
  } else {
    df <- as.data.frame(mat) |> tibble::rownames_to_column("row") |>
      tidyr::pivot_longer(-row, names_to="col", values_to="val")
    p <- ggplot(df, aes(col, row, fill=val)) + geom_tile() + theme_bw() +
      theme(axis.text.x=element_text(angle=60, hjust=1)) + labs(x=NULL,y=NULL,title=title)
    if (MAKE_PDF) ggsave(out_pdf, p, width=10, height=max(4,0.12*nrow(mat)+2))
    if (MAKE_PNG) ggsave(out_png, p, width=10, height=max(4,0.12*nrow(mat)+2), dpi=300)
  }
}

plot_groupwise_scatter <- function(df, title, xlab, ylab, out_pdf, out_png) {
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

plot_overall_scatter <- function(df, title, xlab, ylab, out_pdf, out_png) {
  df <- df |> filter(is.finite(x_used), is.finite(y_used))
  st <- spearman_fast(df$x_used, df$y_used)
  lab <- sprintf("Spearman rho=%.3f\np=%.2g\nn=%d", st["rho"], st["p"], as.integer(st["n"]))

  p <- ggplot(df, aes(x=x_used, y=y_used)) +
    geom_point(aes(color=group), size=2.6, alpha=0.9) +
    geom_smooth(aes(group=1), method="lm", se=TRUE, color="black", fill="grey70", alpha=0.25) +
    scale_color_manual(values=COL_MAP, breaks=names(LAB_MAP), labels=LAB_MAP) +
    theme_bw() +
    labs(title=title, x=xlab, y=ylab, color="group") +
    annotate("text", x=Inf, y=Inf, label=lab, hjust=1.05, vjust=1.2, size=4)

  if (MAKE_PDF) ggsave(out_pdf, p, width=7.8, height=6.0)
  if (MAKE_PNG) ggsave(out_png, p, width=7.8, height=6.0, dpi=300)
}

## -------------------------
## 5) 16S: Build Enterococcus species matrix
## -------------------------
dir_init(OUT_ROOT)
write_tsv0(meta_std, file.path(OUT_ROOT, "00_manifest", "metadata_model_40.tsv"))
write_tsv0(tibble(INPUT_16S=INPUT_16S, TPM_FILE=TPM_FILE, CHERRY_FILE=CHERRY_FILE),
           file.path(OUT_ROOT, "00_manifest", "inputs.tsv"))

tab16 <- fread_smart(INPUT_16S, sep="\t")

# sample cols intersection
sample_cols_16 <- names(tab16)[sapply(names(tab16), is_sample_col)]
sample_cols_16 <- intersect(sample_cols_16, sample_order)
if (length(sample_cols_16) == 0) stop("16S table: no overlapping PDxx sample columns with metadata.")

# Identify taxonomy/genus/species
col_tax <- if ("taxonomy" %in% names(tab16)) "taxonomy" else NA_character_
col_genus <- names(tab16)[str_to_lower(names(tab16)) %in% c("genus","g","genus_name")][1]
col_species <- names(tab16)[str_to_lower(names(tab16)) %in% c("species","s","species_name")][1]

if (!is.na(col_genus) && !is.na(col_species)) {
  tab16$genus <- as.character(tab16[[col_genus]])
  tab16$species <- as.character(tab16[[col_species]])
} else if (!is.na(col_tax)) {
  tab16$genus <- vapply(tab16[[col_tax]], parse_tax_rank, character(1), rank_char="g")
  tab16$species <- vapply(tab16[[col_tax]], parse_tax_rank, character(1), rank_char="s")
} else {
  stop("16S table must contain either (Genus & Species columns) OR a 'taxonomy' column.")
}

s16_ent <- tab16 |>
  mutate(genus = str_trim(genus),
         species = str_trim(species)) |>
  filter(tolower(genus) == "enterococcus") |>
  mutate(species_label = vapply(species, function(s) normalize_species_label("Enterococcus", s), character(1))) |>
  select(species_label, all_of(sample_cols_16))

for (cc in sample_cols_16) s16_ent[[cc]] <- suppressWarnings(as.numeric(s16_ent[[cc]]))
s16_ent[is.na(s16_ent)] <- 0

s16_species <- s16_ent |>
  group_by(species_label) |>
  summarise(across(all_of(sample_cols_16), sum), .groups="drop")

mat16 <- s16_species |>
  tibble::column_to_rownames("species_label") |>
  as.matrix()
storage.mode(mat16) <- "numeric"

# filter
prev16 <- rowMeans(mat16 > 0)
tot16 <- rowSums(mat16)
keep16 <- (prev16 >= MIN_PREV_ALL) & (tot16 >= MIN_TOTAL_16S)

diag16 <- tibble(
  feature = rownames(mat16),
  prev_all = prev16,
  total_count = tot16,
  keep = keep16
) |> arrange(desc(keep), desc(prev_all), desc(total_count))

write_tsv0(diag16, file.path(OUT_ROOT, "tables", "16S_EnterococcusSpecies_filter_diagnostics.tsv"))
write_tsv0(as.data.frame(mat16) |> tibble::rownames_to_column("feature"),
           file.path(OUT_ROOT, "tables", "16S_EnterococcusSpecies_counts.tsv"))

mat16_f <- mat16[keep16, , drop=FALSE]
# CLR on filtered
clr16 <- clr_transform(mat16_f, pseudocount=PSEUDOCOUNT_16S)
clr16 <- clr16[, sample_order, drop=FALSE]
write_tsv0(as.data.frame(clr16) |> tibble::rownames_to_column("feature"),
           file.path(OUT_ROOT, "tables", "16S_EnterococcusSpecies_CLR.tsv"))

## -------------------------
## 6) Virome: Build Host_species matrices (All + Enterococcus subset)
## -------------------------
tpm <- fread_smart(TPM_FILE, sep="\t")
if (!("id" %in% names(tpm))) stop("vOTU_abundance_matrix_TPM.tsv must contain column 'id'.")

sample_cols_v <- setdiff(names(tpm), "id")
sample_cols_v <- intersect(sample_cols_v, sample_order)
if (length(sample_cols_v) == 0) stop("Virome TPM table: no overlapping sample columns with metadata.")

cherry <- fread_smart(CHERRY_FILE, sep="\t")
if (!all(c("Accession","Host") %in% names(cherry))) stop("cherry_prediction.tsv must contain columns: Accession, Host")

v <- tpm |>
  left_join(cherry, by=c("id"="Accession"))

parse_host <- function(host_str) {
  if (is.na(host_str)) return(list(genus=NA_character_, host_species_label=NA_character_))
  s <- str_trim(as.character(host_str))
  m <- str_match(s, "^(species|genus):\\s*(.+)$")
  if (is.na(m[1,1])) return(list(genus=NA_character_, host_species_label=NA_character_))
  rank <- m[1,2]; taxon <- str_trim(m[1,3])
  genus <- str_split(taxon, "\\s+", simplify=TRUE)[1]
  host_species_label <- if (rank == "species") taxon else paste0(genus, " [genus]")
  list(genus=genus, host_species_label=host_species_label)
}
tmp <- lapply(v$Host, parse_host)
v$host_genus <- vapply(tmp, `[[`, character(1), "genus")
v$host_species_label <- vapply(tmp, `[[`, character(1), "host_species_label")

# (Y-all) all host_species labels
v_all <- v |>
  filter(!is.na(host_species_label) & host_species_label != "") |>
  select(host_species_label, all_of(sample_cols_v))
for (cc in sample_cols_v) v_all[[cc]] <- suppressWarnings(as.numeric(v_all[[cc]]))
v_all[is.na(v_all)] <- 0

hs_all <- v_all |>
  group_by(host_species_label) |>
  summarise(across(all_of(sample_cols_v), sum), .groups="drop") |>
  mutate(feature = paste0(host_species_label, "_virus")) |>
  select(feature, all_of(sample_cols_v))

matV_all <- hs_all |> tibble::column_to_rownames("feature") |> as.matrix()
storage.mode(matV_all) <- "numeric"

prevV_all <- rowMeans(matV_all > 0)
totV_all <- rowSums(matV_all)
keepV_all <- (prevV_all >= MIN_PREV_ALL) & (totV_all >= MIN_TOTAL_VIROME)

diagV_all <- tibble(
  feature = rownames(matV_all),
  prev_all = prevV_all,
  total_TPM = totV_all,
  keep = keepV_all
) |> arrange(desc(keep), desc(prev_all), desc(total_TPM))
write_tsv0(diagV_all, file.path(OUT_ROOT, "tables", "Virome_HostSpecies_ALL_filter_diagnostics.tsv"))
write_tsv0(as.data.frame(matV_all) |> tibble::rownames_to_column("feature"),
           file.path(OUT_ROOT, "tables", "Virome_HostSpecies_ALL_TPM.tsv"))

matV_all_f <- matV_all[keepV_all, , drop=FALSE]
clrV_all <- clr_transform(matV_all_f, pseudocount=PSEUDOCOUNT_VIROME)
clrV_all <- clrV_all[, sample_order, drop=FALSE]
write_tsv0(as.data.frame(clrV_all) |> tibble::rownames_to_column("feature"),
           file.path(OUT_ROOT, "tables", "Virome_HostSpecies_ALL_CLR.tsv"))

# (Y-ent) Enterococcus-host subset
v_ent <- v |>
  filter(host_genus == "Enterococcus" & !is.na(host_species_label) & host_species_label != "") |>
  select(host_species_label, all_of(sample_cols_v))
for (cc in sample_cols_v) v_ent[[cc]] <- suppressWarnings(as.numeric(v_ent[[cc]]))
v_ent[is.na(v_ent)] <- 0

hs_ent <- v_ent |>
  group_by(host_species_label) |>
  summarise(across(all_of(sample_cols_v), sum), .groups="drop") |>
  mutate(feature = paste0(host_species_label, "_virus")) |>
  select(feature, all_of(sample_cols_v))

matV_ent <- hs_ent |> tibble::column_to_rownames("feature") |> as.matrix()
storage.mode(matV_ent) <- "numeric"

prevV_ent <- rowMeans(matV_ent > 0)
totV_ent <- rowSums(matV_ent)
keepV_ent <- (prevV_ent >= MIN_PREV_ALL) & (totV_ent >= MIN_TOTAL_VIROME)

diagV_ent <- tibble(
  feature = rownames(matV_ent),
  prev_all = prevV_ent,
  total_TPM = totV_ent,
  keep = keepV_ent
) |> arrange(desc(keep), desc(prev_all), desc(total_TPM))
write_tsv0(diagV_ent, file.path(OUT_ROOT, "tables", "Virome_HostSpecies_Enterococcus_filter_diagnostics.tsv"))
write_tsv0(as.data.frame(matV_ent) |> tibble::rownames_to_column("feature"),
           file.path(OUT_ROOT, "tables", "Virome_HostSpecies_Enterococcus_TPM.tsv"))

matV_ent_f <- matV_ent[keepV_ent, , drop=FALSE]
clrV_ent <- clr_transform(matV_ent_f, pseudocount=PSEUDOCOUNT_VIROME)
clrV_ent <- clrV_ent[, sample_order, drop=FALSE]
write_tsv0(as.data.frame(clrV_ent) |> tibble::rownames_to_column("feature"),
           file.path(OUT_ROOT, "tables", "Virome_HostSpecies_Enterococcus_CLR.tsv"))

## -------------------------
## 7) Correlation runner: X=16S Enterococcus species, Y=Host_species_virus
## -------------------------
run_pair_corr <- function(out_root, X_clr, Y_clr, pair_name, heatmap_mode=c("minp_cols_p05","all_cols")) {
  heatmap_mode <- match.arg(heatmap_mode)
  dir_init(out_root)

  # save manifests
  write_tsv0(meta_std, file.path(out_root, "00_manifest", "metadata_model_40.tsv"))
  write_tsv0(tibble(pair=pair_name, heatmap_mode=heatmap_mode),
             file.path(out_root, "00_manifest", "pair_info.tsv"))

  for (m in MODELS) {
    out_m <- file.path(out_root, m)
    dir_init(out_m)

    # complete-case samples for model covariates
    if (m == "unadjusted") {
      keep_idx <- rep(TRUE, nrow(meta_std))
    } else if (m == "primary") {
      keep_idx <- complete.cases(meta_std[, c("group","age","bmi","disease_duration","constipation")])
    } else {
      keep_idx <- complete.cases(meta_std[, c("group","age","bmi","disease_duration","constipation","levodopa_dose")])
    }
    keep_samples <- meta_std$sample_id[keep_idx]
    n_used <- length(keep_samples)
    if (n_used < 20) {
      writeLines(paste0("Too few complete-case samples for model ", m, ": n=", n_used),
                 con=file.path(out_m, "logs", "note.txt"))
      next
    }

    X_use <- X_clr[, keep_samples, drop=FALSE]
    Y_use <- Y_clr[, keep_samples, drop=FALSE]

    # residualize for adjusted models
    if (m != "unadjusted") {
      Xdes <- build_design(meta_std[keep_idx, , drop=FALSE], m)
      # residualize samples x features
      Xres <- t(residualize_matrix_safe(t(X_use), Xdes))
      Yres <- t(residualize_matrix_safe(t(Y_use), Xdes))

      rownames(Xres) <- rownames(X_use); colnames(Xres) <- colnames(X_use)
      rownames(Yres) <- rownames(Y_use); colnames(Yres) <- colnames(Y_use)

      X_use <- Xres
      Y_use <- Yres
    }

    # save per-sample matrices used (reproducible)
    write_tsv0(as.data.frame(X_use) |> tibble::rownames_to_column("feature"),
               file.path(out_m, "per_sample", "X_values_used_features_by_sample.tsv"))
    write_tsv0(as.data.frame(Y_use) |> tibble::rownames_to_column("feature"),
               file.path(out_m, "per_sample", "Y_values_used_features_by_sample.tsv"))
    write_tsv0(meta_std[keep_idx, ], file.path(out_m, "per_sample", "sample_annotation_used.tsv"))

    # correlation matrices (loop: sizes are manageable for Enterococcus species)
    X_feats <- rownames(X_use); Y_feats <- rownames(Y_use)
    rho <- matrix(NA_real_, nrow=length(X_feats), ncol=length(Y_feats), dimnames=list(X_feats, Y_feats))
    pmat <- rho
    nmat <- rho

    for (i in seq_along(X_feats)) {
      xi <- as.numeric(X_use[X_feats[i], ])
      for (j in seq_along(Y_feats)) {
        yj <- as.numeric(Y_use[Y_feats[j], ])
        st <- spearman_fast(xi, yj)
        rho[i, j] <- as.numeric(st["rho"])
        pmat[i, j] <- as.numeric(st["p"])
        nmat[i, j] <- as.integer(st["n"])
      }
    }

    # global BH across all pairs
    p_vec <- as.vector(pmat)
    q_vec <- p.adjust(p_vec, method="BH")
    qmat <- matrix(q_vec, nrow=nrow(pmat), ncol=ncol(pmat), dimnames=dimnames(pmat))

    # long table
    long <- as.data.frame(as.table(rho), stringsAsFactors=FALSE)
    colnames(long) <- c("X_feature","Y_feature","rho")
    long$p_value <- as.vector(pmat)
    long$q_value <- as.vector(qmat)
    long$n_used <- as.vector(nmat)

    write_tsv0(long |> arrange(q_value),
               file.path(out_m, "tables", "pairwise_correlations_long.tsv"))
    write_tsv0(as.data.frame(rho) |> tibble::rownames_to_column("X_feature"),
               file.path(out_m, "tables", "correlation_matrix_rho.tsv"))
    write_tsv0(as.data.frame(pmat) |> tibble::rownames_to_column("X_feature"),
               file.path(out_m, "tables", "correlation_matrix_p.tsv"))
    write_tsv0(as.data.frame(qmat) |> tibble::rownames_to_column("X_feature"),
               file.path(out_m, "tables", "correlation_matrix_q.tsv"))
    write_tsv0(as.data.frame(nmat) |> tibble::rownames_to_column("X_feature"),
               file.path(out_m, "tables", "correlation_matrix_n.tsv"))

    # heatmap subset selection
    col_minp <- apply(pmat, 2, function(x) suppressWarnings(min(x, na.rm=TRUE)))
    col_minp_tbl <- tibble(Y_feature=names(col_minp), min_p=as.numeric(col_minp)) |>
      arrange(min_p)
    write_tsv0(col_minp_tbl, file.path(out_m, "tables", "heatmap_column_minp.tsv"))

    if (heatmap_mode == "minp_cols_p05") {
      sel_cols <- col_minp_tbl |> filter(is.finite(min_p) & min_p < P_CUTOFF) |> pull(Y_feature)
      if (length(sel_cols) == 0) sel_cols <- col_minp_tbl |> head(HEATMAP_MAX_COLS) |> pull(Y_feature)
      if (length(sel_cols) > HEATMAP_MAX_COLS) sel_cols <- col_minp_tbl |> filter(Y_feature %in% sel_cols) |> head(HEATMAP_MAX_COLS) |> pull(Y_feature)
      sel_cols <- intersect(sel_cols, colnames(rho))
      write_tsv0(tibble(Y_feature=sel_cols), file.path(out_m, "tables", "heatmap_selected_columns.tsv"))

      sub_mat <- rho[, sel_cols, drop=FALSE]  # keep ALL rows
    } else {
      # all columns (after prevalence filter)
      sub_mat <- rho
      write_tsv0(tibble(Y_feature=colnames(rho)), file.path(out_m, "tables", "heatmap_selected_columns.tsv"))
    }

    write_tsv0(as.data.frame(sub_mat) |> tibble::rownames_to_column("X_feature"),
               file.path(out_m, "tables", "heatmap_subset_rho.tsv"))

    plot_heatmap(sub_mat,
                 out_pdf=file.path(out_m, "figures_pdf", "heatmap_rho.pdf"),
                 out_png=file.path(out_m, "figures_png", "heatmap_rho.png"),
                 title=paste0(pair_name, " (", m, "), n=", n_used))

    writeLines(paste0("Done: ", pair_name, " | ", m, " | n_used=", n_used,
                      " | X=", nrow(X_use), " features, Y=", nrow(Y_use), " features"),
               con=file.path(out_m, "logs", "done.txt"))
  }
}

## -------------------------
## 8) Run two correlation modules
## -------------------------
OUT_ALL <- file.path(OUT_ROOT, "03_Corr_EntSpecies_vs_AllHostSpecies")
OUT_ENT <- file.path(OUT_ROOT, "04_Corr_EntSpecies_vs_EntHostSpecies")

run_pair_corr(OUT_ALL, clr16, clrV_all,
              pair_name="16S_EnterococcusSpecies__Virome_HostSpecies_ALL",
              heatmap_mode="minp_cols_p05")

run_pair_corr(OUT_ENT, clr16, clrV_ent,
              pair_name="16S_EnterococcusSpecies__Virome_HostSpecies_Enterococcus",
              heatmap_mode="all_cols")

## -------------------------
## 9) Scatter plots (matched species pairs), same style as old script
## -------------------------
if (DO_SCATTER) {
  OUT_SC <- file.path(OUT_ROOT, "05_Scatter_MatchedSpeciesPairs")
  dir_init(OUT_SC)

  # match by removing "_virus" suffix; case-insensitive exact match
  x_names <- rownames(clr16)
  y_names <- rownames(clrV_all)  # use ALL host species for matching
  y_base <- str_replace(y_names, "_virus$", "")
  # standardize spaces
  x_key <- str_squish(x_names)
  y_key <- str_squish(y_base)

  common <- intersect(tolower(x_key), tolower(y_key))
  if (length(common) == 0) {
    writeLines("No matched species names between 16S species and virome host_species labels.", con=file.path(OUT_SC, "logs", "note.txt"))
  } else {
    # build mapping
    map <- tibble(
      x_feature = x_names[match(common, tolower(x_key))],
      y_feature = y_names[match(common, tolower(y_key))]
    )
    write_tsv0(map, file.path(OUT_SC, "00_manifest", "matched_pairs.tsv"))

    get_xy_by_model <- function(x_vec, y_vec, model_name, meta_df) {
      if (model_name == "unadjusted") return(list(x=x_vec, y=y_vec, keep=rep(TRUE, nrow(meta_df))))
      if (model_name == "primary") {
        keep <- complete.cases(meta_df[, c("group","age","bmi","disease_duration","constipation")])
        Xdes <- build_design(meta_df[keep, , drop=FALSE], "primary")
      } else {
        keep <- complete.cases(meta_df[, c("group","age","bmi","disease_duration","constipation","levodopa_dose")])
        Xdes <- build_design(meta_df[keep, , drop=FALSE], "sensitivity_add_dose")
      }
      xr <- rep(NA_real_, length(x_vec)); yr <- rep(NA_real_, length(y_vec))
      fitx <- lm.fit(Xdes, x_vec[keep]); fity <- lm.fit(Xdes, y_vec[keep])
      xr[keep] <- fitx$residuals
      yr[keep] <- fity$residuals
      list(x=xr, y=yr, keep=keep)
    }

    for (m in MODELS) {
      out_m <- file.path(OUT_SC, m)
      dir_init(out_m)

      for (k in seq_len(nrow(map))) {
        sp <- map$x_feature[k]
        vp <- map$y_feature[k]

        x_raw <- as.numeric(clr16[sp, sample_order])
        y_raw <- as.numeric(clrV_all[vp, sample_order])

        xy <- get_xy_by_model(x_raw, y_raw, m, meta_std)

        df <- tibble(
          sample_id = sample_order,
          group = meta_std$group,
          x_used = xy$x,
          y_used = xy$y
        )

        # save per-sample
        write_tsv0(df, file.path(out_m, "per_sample", paste0("scatter_data_", sp, "__", vp, ".tsv")))

        # overall
        plot_overall_scatter(
          df,
          title=paste0("Overall | ", m, " | X=16S ", sp, " ; Y=Phage host ", str_replace(vp, "_virus$", "")),
          xlab=if (m=="unadjusted") paste0("16S ", sp, " (CLR)") else paste0("16S ", sp, " (resid, ", m, ")"),
          ylab=if (m=="unadjusted") paste0("Phage host ", str_replace(vp, "_virus$", ""), " (CLR)") else paste0("Phage host ", str_replace(vp, "_virus$", ""), " (resid, ", m, ")"),
          out_pdf=file.path(out_m, "figures_pdf", paste0("scatter_overall__", sp, "__", vp, ".pdf")),
          out_png=file.path(out_m, "figures_png", paste0("scatter_overall__", sp, "__", vp, ".png"))
        )

        # groupwise
        plot_groupwise_scatter(
          df,
          title=paste0("Groupwise | ", m, " | X=16S ", sp, " ; Y=Phage host ", str_replace(vp, "_virus$", "")),
          xlab=if (m=="unadjusted") paste0("16S ", sp, " (CLR)") else paste0("16S ", sp, " (resid, ", m, ")"),
          ylab=if (m=="unadjusted") paste0("Phage host ", str_replace(vp, "_virus$", ""), " (CLR)") else paste0("Phage host ", str_replace(vp, "_virus$", ""), " (resid, ", m, ")"),
          out_pdf=file.path(out_m, "figures_pdf", paste0("scatter_groupwise__", sp, "__", vp, ".pdf")),
          out_png=file.path(out_m, "figures_png", paste0("scatter_groupwise__", sp, "__", vp, ".png"))
        )
      }
    }
  }
}

sink(file.path(OUT_ROOT, "00_manifest", "sessionInfo.txt"))
print(sessionInfo())
sink()

cat("All done -> ", OUT_ROOT, "\n")