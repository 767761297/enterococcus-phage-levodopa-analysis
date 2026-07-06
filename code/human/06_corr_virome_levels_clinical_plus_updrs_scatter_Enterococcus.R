## =========================
## 11_corr_virome_levels_clinical_plus_updrs_scatter_Enterococcus.R
## Virome correlation full suite + UPDRS scatterplots for Enterococcus phage
##
## Inputs in working dir:
##   metadata.csv
##   vOTU_abundance_matrix_TPM.tsv
##   Virus_taxonomy.tsv          (2 cols, no header)
##   cherry_prediction.tsv
##
## Outputs:
##   results_virome_correlation/
##     Host/ Host_species/ ViralTax_phylum ... ViralTax_species/
##     Each level has unadjusted/primary/sensitivity_add_dose subfolders
##     Tables: long + matrices + per-sample reproducible data
##     Figures: heatmaps (PDF) + UPDRS scatter (Host levels only, PDF one figure per file)
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

## ---------- Inputs ----------
META_FILE <- "metadata.csv"
TPM_FILE  <- "vOTU_abundance_matrix_TPM.tsv"
TAX_FILE  <- "Virus_taxonomy.tsv"        # 2 cols, no header
HOST_FILE <- "cherry_prediction.tsv"

## ---------- Output ----------
OUT_ROOT <- "results_virome_correlation"
dir.create(OUT_ROOT, showWarnings = FALSE)

## ---------- Global params ----------
GROUP_LEVELS <- c("PD_GR","PD_MR")
VIRAL_RANKS <- c("phylum","class","order","family","genus","species")

# Filtering (independent)
MIN_PREVALENCE_ALL <- 0.10
MIN_TOTAL_TPM <- 0
MIN_PREV_BY_GROUP <- NA_real_  # optional e.g. 0.05
MIN_MEAN_RA <- NA_real_        # optional e.g. 1e-5

# CLR
PSEUDOCOUNT <- 0.01

# Heatmap selection
P_CUTOFF <- 0.05
FALLBACK_TOP_N <- 50
MAKE_PDF <- TRUE
MAKE_PNG <- TRUE
CLUSTER_ROWS <- TRUE
CLUSTER_COLS <- FALSE

# Scatterplots (Host only)
MAKE_UPDRS_SCATTER <- TRUE
COL_MAP <- c(PD_GR = "#7CAE00", PD_MR = "#F8766D")
LAB_MAP <- c(PD_GR = "GR", PD_MR = "MR")

MODELS <- c("unadjusted","primary","sensitivity_add_dose")

## ---------- Helpers ----------
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
  dir.create(file.path(path, "00_manifest"), showWarnings = FALSE)
  dir.create(file.path(path, "per_sample"), showWarnings = FALSE)
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
safe_mode <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != "" & x != "NA"]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
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

# CLR: features x samples
clr_transform <- function(mat_feat_by_samp, pseudocount = 0.01) {
  x <- mat_feat_by_samp + pseudocount
  x <- sweep(x, 2, colSums(x), "/")
  lx <- log(x)
  sweep(lx, 2, colMeans(lx), "-")
}
mean_ra_fast <- function(mat_feat_by_samp) {
  tot <- colSums(mat_feat_by_samp)
  w <- 1 / tot
  as.numeric(mat_feat_by_samp %*% w) / length(w)
}
filter_matrix <- function(mat_feat_by_samp, meta,
                          min_prev_all, min_total,
                          min_prev_by_group = NA_real_,
                          min_mean_ra = NA_real_) {
  prev_all <- rowMeans(mat_feat_by_samp > 0)
  tot <- rowSums(mat_feat_by_samp)
  keep <- (prev_all >= min_prev_all) & (tot >= min_total)

  if (!is.na(min_prev_by_group)) {
    idx_gr <- which(meta$group=="PD_GR")
    idx_mr <- which(meta$group=="PD_MR")
    prev_gr <- rowMeans(mat_feat_by_samp[, idx_gr, drop=FALSE] > 0)
    prev_mr <- rowMeans(mat_feat_by_samp[, idx_mr, drop=FALSE] > 0)
    keep <- keep & (prev_gr >= min_prev_by_group) & (prev_mr >= min_prev_by_group)
  } else {
    prev_gr <- prev_mr <- rep(NA_real_, nrow(mat_feat_by_samp))
  }

  if (!is.na(min_mean_ra)) {
    mra <- mean_ra_fast(mat_feat_by_samp)
    keep <- keep & (mra >= min_mean_ra)
  } else {
    mra <- rep(NA_real_, nrow(mat_feat_by_samp))
  }

  list(
    mat = mat_feat_by_samp[keep, , drop=FALSE],
    diag = tibble(feature=rownames(mat_feat_by_samp),
                  prev_all=prev_all, prev_gr=prev_gr, prev_mr=prev_mr,
                  total=tot, mean_ra=mra, keep=keep),
    summary = tibble(n_total=nrow(mat_feat_by_samp), n_keep=sum(keep),
                     min_prev_all=min_prev_all, min_total=min_total,
                     min_prev_by_group=min_prev_by_group, min_mean_ra=min_mean_ra)
  )
}
aggregate_sum <- function(mat, map_named) {
  grp <- map_named[rownames(mat)]
  rowsum(mat, group = grp, reorder = TRUE)
}

# Virus taxonomy parser
parse_virus_tax <- function(x) {
  parts <- str_split(as.character(x), ";", simplify = TRUE)
  get_rank <- function(row, prefix) {
    y <- row[str_detect(row, paste0("^\\s*", prefix))]
    if (length(y)==0) return(NA_character_)
    str_replace(str_trim(y[1]), paste0("^", prefix), "")
  }
  tibble(
    realm   = apply(parts, 1, get_rank, prefix="r__"),
    kingdom = apply(parts, 1, get_rank, prefix="k__"),
    phylum  = apply(parts, 1, get_rank, prefix="p__"),
    class   = apply(parts, 1, get_rank, prefix="c__"),
    order   = apply(parts, 1, get_rank, prefix="o__"),
    family  = apply(parts, 1, get_rank, prefix="f__"),
    genus   = apply(parts, 1, get_rank, prefix="g__"),
    species = apply(parts, 1, get_rank, prefix="s__")
  )
}

# Host parser: outputs host_genus, host_species
parse_host_fields <- function(host_raw) {
  host_raw <- as.character(host_raw)
  host_raw[host_raw %in% c("-", "", "NA")] <- NA_character_

  rank <- ifelse(!is.na(host_raw) & str_detect(host_raw, ":"), str_split_fixed(host_raw, ":", 2)[,1], NA_character_)
  name <- ifelse(!is.na(host_raw) & str_detect(host_raw, ":"), str_split_fixed(host_raw, ":", 2)[,2], host_raw)
  name <- str_trim(name)

  host_species <- ifelse(!is.na(rank) & rank=="species" & !is.na(name) & name!="", name, NA_character_)
  host_genus <- ifelse(!is.na(rank) & rank=="genus" & !is.na(name) & name!="", name,
                       ifelse(!is.na(host_species) & str_detect(host_species, "\\s+"),
                              str_split_fixed(host_species, "\\s+", 2)[,1],
                              NA_character_))
  host_genus <- ifelse(!is.na(host_genus) & host_genus!="", host_genus, "Unassigned_host")

  tibble(host_rank=rank, host_name=name, host_genus=host_genus, host_species=host_species)
}

# Residualization (matrix + vector)
residualize_matrix <- function(Y_samp_by_feat, X_design) {
  qrX <- qr(X_design)
  qr.resid(qrX, Y_samp_by_feat)
}
residualize_vector <- function(y, X_design) {
  fit <- lm.fit(X_design, y)
  fit$residuals
}

# Spearman and Wilcoxon (rank-biserial)
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
  z <- c(x1, x0)
  r <- rank(z, ties.method = "average")
  Wrank <- sum(r[seq_len(n1)])
  U <- Wrank - n1*(n1+1)/2
  rb <- (2*U)/(n1*n0) - 1
  rb <- max(-1, min(1, rb))
  c(rb=rb, p=p, n=n1+n0)
}

plot_heatmap_effect <- function(effect_mat, out_pdf, out_png) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    df <- as.data.frame(effect_mat) |> rownames_to_column("feature") |>
      pivot_longer(-feature, names_to="covariate", values_to="effect")
    df$feature <- factor(df$feature, levels = rev(rownames(effect_mat)))
    p <- ggplot(df, aes(covariate, feature, fill=effect)) + geom_tile() + theme_bw() + labs(x=NULL,y=NULL)
    if (MAKE_PDF) ggsave(out_pdf, p, width=9, height=max(4, 0.18*nrow(effect_mat)+2))
    if (MAKE_PNG) ggsave(out_png, p, width=9, height=max(4, 0.18*nrow(effect_mat)+2), dpi=300)
    return(invisible(NULL))
  }
  cr <- CLUSTER_ROWS && nrow(effect_mat) >= 2
  cc <- CLUSTER_COLS && ncol(effect_mat) >= 2
  if (MAKE_PDF) {
    pdf(out_pdf, width=9, height=max(4, 0.18*nrow(effect_mat)+2))
    pheatmap::pheatmap(effect_mat, cluster_rows=cr, cluster_cols=cc, show_colnames=TRUE, fontsize_row=8)
    dev.off()
  }
  if (MAKE_PNG) {
    png(out_png, width=2200, height=max(1400, 40*nrow(effect_mat)+600), res=300)
    pheatmap::pheatmap(effect_mat, cluster_rows=cr, cluster_cols=cc, show_colnames=TRUE, fontsize_row=8)
    dev.off()
  }
}

# Scatter residualization for UPDRS: within-group or overall
resid_on_covars <- function(df, value, covars) {
  keep <- complete.cases(df[, covars, drop=FALSE]) & is.finite(value)
  out <- rep(NA_real_, length(value))
  if (sum(keep) < 6) return(out)
  X <- model.matrix(as.formula(paste("~", paste(covars, collapse=" + "))), data=df[keep,,drop=FALSE])
  fit <- lm.fit(X, value[keep])
  out[keep] <- fit$residuals
  out
}
make_xy_scatter <- function(df, x_raw, y_raw, model, mode=c("overall","within_group")) {
  mode <- match.arg(mode)
  cov_primary <- c("age","bmi","disease_duration","constipation")
  cov_sens <- c(cov_primary, "levodopa_dose")

  if (model=="unadjusted") {
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
      df$x_used[idx] <- resid_on_covars(df[idx,,drop=FALSE], x_raw[idx], covars)
      df$y_used[idx] <- resid_on_covars(df[idx,,drop=FALSE], y_raw[idx], covars)
    }
  }
  df
}
plot_groupwise_scatter <- function(df, title, xlab, ylab, out_pdf, out_png) {
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

## ---------- Load metadata ----------
meta_raw <- read_csv_smart(META_FILE) |> janitor::clean_names(replace=c("\u03bc"="u"))
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
    group = factor(group, levels=GROUP_LEVELS),
    age = suppressWarnings(as.numeric(.data[[age_col]])),
    bmi = suppressWarnings(as.numeric(.data[[bmi_col]])),
    disease_duration = suppressWarnings(as.numeric(.data[[dur_col]])),
    constipation = factor(as.integer(.data[[con_col]]), levels=c(0,1), labels=c("No","Yes")),
    levodopa_dose = suppressWarnings(as.numeric(.data[[dose_col]])),
    updrs3_reduction_pct = parse_percent(.data[[updrs_col]])
  ) |>
  arrange(group, sample_id)

stopifnot(nrow(meta)==40)
sample_order <- meta$sample_id

## Covariate table
COVARS <- tibble::tribble(
  ~covar, ~type,
  "age", "numeric",
  "bmi", "numeric",
  "disease_duration", "numeric",
  "constipation", "binary",
  "levodopa_dose", "numeric",
  "updrs3_reduction_pct", "numeric"
)

## ---------- Read TPM matrix ----------
tpm_dt <- fread(TPM_FILE, sep="\t", header=TRUE, data.table=FALSE, check.names=FALSE)
stopifnot("id" %in% names(tpm_dt))

sample_cols <- setdiff(names(tpm_dt), "id")
sample_cols <- sample_cols[str_detect(sample_cols, "^PD\\d+$")]
keep_samples <- intersect(sample_order, sample_cols)
stopifnot(length(keep_samples)==40)
keep_samples <- sample_order

mat_votu <- tpm_dt[, keep_samples, drop=FALSE]
mat_votu[] <- lapply(mat_votu, as.numeric)
mat_votu <- as.matrix(mat_votu)
rownames(mat_votu) <- as.character(tpm_dt$id)

## ---------- Read annotations ----------
tax_dt <- fread(TAX_FILE, sep="\t", header=FALSE, data.table=FALSE)
colnames(tax_dt) <- c("votu_id","taxonomy_raw")
tax_tbl <- tax_dt |>
  mutate(votu_id = as.character(votu_id)) |>
  bind_cols(parse_virus_tax(tax_dt$taxonomy_raw)) |>
  mutate(across(all_of(VIRAL_RANKS), ~fill_unknown(.x, cur_column()))) |>
  mutate(realm = fill_unknown(realm, "realm"),
         kingdom = fill_unknown(kingdom, "kingdom"))

host_dt <- fread(HOST_FILE, sep="\t", header=TRUE, data.table=FALSE)
# robust column access by name; most likely as in your files
host_tbl <- host_dt |>
  transmute(
    votu_id = as.character(Accession),
    Host_raw = as.character(Host),
    cherry_score = suppressWarnings(as.numeric(CHERRYScore)),
    method = as.character(Method),
    host_ncbi_lineage_raw = as.character(Host_NCBI_lineage),
    host_gtdb_lineage_raw = as.character(Host_GTDB_lineage)
  )
host_tbl <- bind_cols(host_tbl, parse_host_fields(host_tbl$Host_raw)) |>
  mutate(host_genus = fill_unknown(host_genus, "host"))

anno_votu <- tibble(votu_id = rownames(mat_votu)) |>
  left_join(tax_tbl, by="votu_id") |>
  left_join(host_tbl, by="votu_id") |>
  mutate(
    host_genus = fill_unknown(host_genus, "host"),
    host_species = as.character(host_species),
    method = ifelse(is.na(method) | method=="", "Unspecified", method)
  )

## ---------- Base filter on vOTU (backbone) ----------
f0 <- filter_matrix(mat_votu, meta,
                    min_prev_all = MIN_PREVALENCE_ALL,
                    min_total = MIN_TOTAL_TPM,
                    min_prev_by_group = MIN_PREV_BY_GROUP,
                    min_mean_ra = MIN_MEAN_RA)
mat_votu0 <- f0$mat
anno_votu0 <- anno_votu |> filter(votu_id %in% rownames(mat_votu0))

## ---------- Build level matrices ----------
build_level <- function(level, rank=NULL, mat_votu_in, anno_votu_in) {
  if (level=="Host") {
    map <- setNames(anno_votu_in$host_genus, anno_votu_in$votu_id)
    mat <- aggregate_sum(mat_votu_in, map)
    feat_anno <- anno_votu_in |>
      group_by(host_genus) |>
      summarise(
        n_member_votu = n(),
        member_votu_ids = paste(votu_id, collapse=";"),
        top_method = safe_mode(method),
        mean_cherry_score = if (all(is.na(cherry_score))) NA_real_ else mean(cherry_score, na.rm=TRUE),
        max_cherry_score  = if (all(is.na(cherry_score))) NA_real_ else max(cherry_score, na.rm=TRUE),
        .groups="drop"
      ) |>
      mutate(feature = host_genus,
             level = "Host_genus",
             display_label = host_genus) |>
      select(feature, display_label, level,
             n_member_votu, member_votu_ids, top_method, mean_cherry_score, max_cherry_score)
    return(list(mat=mat, feat_anno=feat_anno))
  }

  if (level=="Host_species") {
    host_species_key <- ifelse(!is.na(anno_votu_in$host_species) & anno_votu_in$host_species!="",
                               anno_votu_in$host_species,
                               paste0(anno_votu_in$host_genus, " [genus]"))
    host_species_key <- fill_unknown(host_species_key, "host_species")
    map <- setNames(host_species_key, anno_votu_in$votu_id)
    mat <- aggregate_sum(mat_votu_in, map)
    feat_anno <- anno_votu_in |>
      mutate(host_species_key = map[votu_id]) |>
      group_by(host_species_key) |>
      summarise(
        n_member_votu = n(),
        member_votu_ids = paste(votu_id, collapse=";"),
        top_host_genus = safe_mode(host_genus),
        top_method = safe_mode(method),
        mean_cherry_score = if (all(is.na(cherry_score))) NA_real_ else mean(cherry_score, na.rm=TRUE),
        max_cherry_score  = if (all(is.na(cherry_score))) NA_real_ else max(cherry_score, na.rm=TRUE),
        .groups="drop"
      ) |>
      mutate(feature = host_species_key,
             level = "Host_species",
             display_label = host_species_key) |>
      select(feature, display_label, level,
             n_member_votu, member_votu_ids, top_host_genus, top_method, mean_cherry_score, max_cherry_score)
    return(list(mat=mat, feat_anno=feat_anno))
  }

  if (level=="ViralTax") {
    stopifnot(!is.null(rank))
    taxon <- fill_unknown(anno_votu_in[[rank]], rank)
    map <- setNames(taxon, anno_votu_in$votu_id)
    mat <- aggregate_sum(mat_votu_in, map)

    feat_anno <- anno_votu_in |>
      mutate(taxon = map[votu_id]) |>
      group_by(taxon) |>
      summarise(
        n_member_votu = n(),
        member_votu_ids = paste(votu_id, collapse=";"),
        realm = safe_mode(realm),
        kingdom = safe_mode(kingdom),
        phylum = safe_mode(phylum),
        class = safe_mode(class),
        order = safe_mode(order),
        family = safe_mode(family),
        genus = safe_mode(genus),
        species = safe_mode(species),
        top_host_genus = safe_mode(host_genus),
        top_host_species = safe_mode(host_species),
        top_method = safe_mode(method),
        mean_cherry_score = if (all(is.na(cherry_score))) NA_real_ else mean(cherry_score, na.rm=TRUE),
        max_cherry_score  = if (all(is.na(cherry_score))) NA_real_ else max(cherry_score, na.rm=TRUE),
        .groups="drop"
      ) |>
      mutate(feature = taxon,
             level = paste0("ViralTax_", rank),
             display_label = taxon) |>
      select(feature, display_label, level, n_member_votu, member_votu_ids,
             realm, kingdom, phylum, class, order, family, genus, species,
             top_host_genus, top_host_species, top_method, mean_cherry_score, max_cherry_score)
    return(list(mat=mat, feat_anno=feat_anno))
  }

  stop("Unknown level")
}

## ---------- Build adjustment design matrices ----------
build_adjust_X_overall <- function(meta, model_name, tested_covar) {
  # overall association adjustment includes group + base covars (exclude tested)
  base <- c("group", "age", "bmi", "disease_duration", "constipation")
  if (model_name == "sensitivity_add_dose") base <- c(base, "levodopa_dose")
  if (model_name == "unadjusted") return(NULL)

  adj <- setdiff(base, tested_covar)
  if (length(adj) == 0) return(model.matrix(~ 1, data = meta))
  model.matrix(as.formula(paste("~", paste(adj, collapse=" + "))), data = meta)
}

build_adjust_X_within_group <- function(meta_g, model_name, tested_covar) {
  base <- c("age", "bmi", "disease_duration", "constipation")
  if (model_name == "sensitivity_add_dose") base <- c(base, "levodopa_dose")
  if (model_name == "unadjusted") return(NULL)

  adj <- setdiff(base, tested_covar)
  if (length(adj) == 0) return(model.matrix(~ 1, data = meta_g))
  model.matrix(as.formula(paste("~", paste(adj, collapse=" + "))), data = meta_g)
}

## ---------- Correlation runner per level ----------
run_level_corr <- function(level_name, level_obj) {
  level_root <- file.path(OUT_ROOT, level_name)
  dir_init(level_root)

  # filter at level
  f <- filter_matrix(level_obj$mat, meta,
                     min_prev_all = MIN_PREVALENCE_ALL,
                     min_total = MIN_TOTAL_TPM,
                     min_prev_by_group = MIN_PREV_BY_GROUP,
                     min_mean_ra = MIN_MEAN_RA)
  mat_f <- f$mat[, sample_order, drop=FALSE]
  feat_anno <- level_obj$feat_anno |> filter(feature %in% rownames(mat_f))

  write_tsv0(f$summary, file.path(level_root, "00_manifest", "filter_summary.tsv"))
  write_tsv0(f$diag, file.path(level_root, "00_manifest", "filter_diagnostics.tsv"))
  write_tsv0(feat_anno, file.path(level_root, "00_manifest", "feature_annotation.tsv"))
  write_tsv0(COVARS, file.path(level_root, "00_manifest", "covariates.tsv"))
  write_tsv0(meta, file.path(level_root, "00_manifest", "metadata_model_40.tsv"))

  # base matrices
  ra_f <- sweep(mat_f, 2, colSums(mat_f), "/")
  clr_f <- clr_transform(mat_f, pseudocount = PSEUDOCOUNT)

  write_tsv0(as.data.frame(mat_f) |> rownames_to_column("feature"),
             file.path(level_root, "per_sample", "TPM_features_by_sample.tsv"))
  write_tsv0(as.data.frame(ra_f) |> rownames_to_column("feature"),
             file.path(level_root, "per_sample", "RA_features_by_sample.tsv"))
  write_tsv0(as.data.frame(clr_f) |> rownames_to_column("feature"),
             file.path(level_root, "per_sample", "CLR_features_by_sample.tsv"))

  for (mname in MODELS) {
    out_dir <- file.path(level_root, mname)
    dir_init(out_dir)
    write_tsv0(tibble(model=mname), file.path(out_dir, "00_manifest", "model.tsv"))

    Y0 <- t(clr_f)  # samples x features
    rownames(Y0) <- colnames(clr_f)
    colnames(Y0) <- rownames(clr_f)

    all_res <- list()
    # matrices: features x covariates
    n_feat <- ncol(Y0)
    eff_mat <- matrix(NA_real_, nrow=n_feat, ncol=nrow(COVARS))
    p_mat   <- matrix(NA_real_, nrow=n_feat, ncol=nrow(COVARS))
    q_mat   <- matrix(NA_real_, nrow=n_feat, ncol=nrow(COVARS))
    rownames(eff_mat) <- colnames(Y0)
    rownames(p_mat) <- colnames(Y0)
    rownames(q_mat) <- colnames(Y0)
    colnames(eff_mat) <- COVARS$covar
    colnames(p_mat) <- COVARS$covar
    colnames(q_mat) <- COVARS$covar

    ## ---- overall association ----
    for (i in seq_len(nrow(COVARS))) {
      cv <- COVARS$covar[i]
      tp <- COVARS$type[i]
      Xadj <- build_adjust_X_overall(meta, mname, cv)

      Y <- Y0
      if (!is.null(Xadj) && mname!="unadjusted") {
        Y <- residualize_matrix(Y0, Xadj)
      }

      cov_raw <- meta[[cv]]
      cov_use <- cov_raw
      if (tp=="numeric" && !is.null(Xadj) && mname!="unadjusted") {
        cov_use <- residualize_vector(cov_raw, Xadj)
      }

      if (tp=="numeric") {
        for (ftr in colnames(Y)) {
          ss <- spearman_fast(Y[, ftr], cov_use)
          all_res[[length(all_res)+1]] <- tibble(
            feature=ftr, covariate=cv, test="spearman",
            effect=ss["rho"], p_value=ss["p"], n=ss["n"]
          )
          eff_mat[ftr, cv] <- as.numeric(ss["rho"])
          p_mat[ftr, cv] <- as.numeric(ss["p"])
        }
      } else {
        g01 <- ifelse(meta[[cv]] == "Yes", 1L, 0L)
        for (ftr in colnames(Y)) {
          ww <- wilcox_rb(Y[, ftr], g01)
          all_res[[length(all_res)+1]] <- tibble(
            feature=ftr, covariate=cv, test="wilcox",
            effect=ww["rb"], p_value=ww["p"], n=ww["n"]
          )
          eff_mat[ftr, cv] <- as.numeric(ww["rb"])
          p_mat[ftr, cv] <- as.numeric(ww["p"])
        }
      }

      # per-sample reproducible inputs used for this covariate/model
      write_tsv0(as.data.frame(t(Y)) |> rownames_to_column("feature"),
                 file.path(out_dir, "per_sample", paste0("feature_values_used_for_", cv, ".tsv")))
      if (tp=="numeric") {
        write_tsv0(tibble(sample_id=meta$sample_id, covariate=cv, value_used=cov_use),
                   file.path(out_dir, "per_sample", paste0("covariate_values_used_for_", cv, ".tsv")))
      }
    }

    res_df <- bind_rows(all_res) |>
      group_by(covariate) |>
      mutate(q_value = p.adjust(p_value, method="BH")) |>
      ungroup() |>
      left_join(feat_anno, by="feature")

    write_tsv0(res_df, file.path(out_dir, "tables", "association_overall_long.tsv"))

    # fill q_mat
    q_tmp <- res_df |> select(feature, covariate, q_value) |>
      pivot_wider(names_from=covariate, values_from=q_value)
    q_df <- as.data.frame(q_tmp); rownames(q_df) <- q_df$feature; q_df$feature <- NULL
    q_mat[rownames(q_df), colnames(q_df)] <- as.matrix(q_df)

    write_tsv0(as.data.frame(eff_mat) |> rownames_to_column("feature"),
               file.path(out_dir, "tables", "effect_matrix.tsv"))
    write_tsv0(as.data.frame(p_mat) |> rownames_to_column("feature"),
               file.path(out_dir, "tables", "p_matrix.tsv"))
    write_tsv0(as.data.frame(q_mat) |> rownames_to_column("feature"),
               file.path(out_dir, "tables", "q_matrix.tsv"))

    # heatmap selection by min p
    min_p <- apply(p_mat, 1, function(x) suppressWarnings(min(x, na.rm=TRUE)))
    sel <- names(min_p)[is.finite(min_p) & min_p < P_CUTOFF]
    if (length(sel)==0) {
      sel <- names(sort(min_p))[seq_len(min(FALLBACK_TOP_N, length(min_p)))]
      writeLines(paste0("No features with min p<", P_CUTOFF, "; fallback top ", length(sel), " by min p."),
                 con=file.path(out_dir, "logs", "heatmap_selection_note.txt"))
    } else {
      sel <- sel[order(min_p[sel])]
      if (length(sel) > FALLBACK_TOP_N) sel <- sel[seq_len(FALLBACK_TOP_N)]
      writeLines(paste0("Selected features: min p<", P_CUTOFF, " (n=", length(sel), "), capped at ", length(sel), "."),
                 con=file.path(out_dir, "logs", "heatmap_selection_note.txt"))
    }
    write_tsv0(tibble(feature=sel, min_p=min_p[sel]),
               file.path(out_dir, "tables", "heatmap_selected_features.tsv"))

    eff_sel <- eff_mat[sel, , drop=FALSE]
    plot_heatmap_effect(eff_sel,
                        out_pdf=file.path(out_dir, "figures_pdf", paste0("heatmap_effect_", level_name, "_", mname, ".pdf")),
                        out_png=file.path(out_dir, "figures_png", paste0("heatmap_effect_", level_name, "_", mname, ".png")))

    ## ---- within-group association ----
    for (gname in GROUP_LEVELS) {
      idx <- which(meta$group == gname)
      meta_g <- meta[idx, , drop=FALSE]
      Yg0 <- t(clr_f[, idx, drop=FALSE])
      rownames(Yg0) <- meta_g$sample_id
      colnames(Yg0) <- rownames(clr_f)

      allg <- list()
      for (i in seq_len(nrow(COVARS))) {
        cv <- COVARS$covar[i]
        tp <- COVARS$type[i]

        Xadj <- build_adjust_X_within_group(meta_g, mname, cv)
        Yg <- Yg0
        if (!is.null(Xadj) && mname!="unadjusted") Yg <- residualize_matrix(Yg0, Xadj)

        cov_raw <- meta_g[[cv]]
        cov_use <- cov_raw
        if (tp=="numeric" && !is.null(Xadj) && mname!="unadjusted") cov_use <- residualize_vector(cov_raw, Xadj)

        if (tp=="numeric") {
          for (ftr in colnames(Yg)) {
            ss <- spearman_fast(Yg[, ftr], cov_use)
            allg[[length(allg)+1]] <- tibble(group=gname, model=mname, feature=ftr, covariate=cv,
                                             test="spearman", effect=ss["rho"], p_value=ss["p"], n=ss["n"])
          }
        } else {
          g01 <- ifelse(meta_g[[cv]]=="Yes", 1L, 0L)
          for (ftr in colnames(Yg)) {
            ww <- wilcox_rb(Yg[, ftr], g01)
            allg[[length(allg)+1]] <- tibble(group=gname, model=mname, feature=ftr, covariate=cv,
                                             test="wilcox", effect=ww["rb"], p_value=ww["p"], n=ww["n"])
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

    ## ---- interaction test (linear model on CLR) ----
    inter_all <- list()
    Yclr <- t(clr_f)  # samples x features
    for (i in seq_len(nrow(COVARS))) {
      cv <- COVARS$covar[i]
      # adjustment covars for model (exclude cv itself)
      base <- c("age","bmi","disease_duration","constipation")
      if (mname=="sensitivity_add_dose") base <- c(base, "levodopa_dose")
      others <- setdiff(base, cv)

      rhs <- paste0("group * ", cv)
      if (length(others) > 0) rhs <- paste(rhs, paste(others, collapse=" + "), sep=" + ")
      fml <- as.formula(paste0("y ~ ", rhs))

      for (ftr in colnames(Yclr)) {
        dfm <- meta |> mutate(y = Yclr[, ftr])
        fit <- try(lm(fml, data=dfm), silent=TRUE)
        if (inherits(fit,"try-error")) next
        coefs <- summary(fit)$coefficients
        iterm <- grep(paste0("^groupPD_MR:", cv, "$|^", cv, ":groupPD_MR$"),
                      rownames(coefs), value=TRUE)
        if (length(iterm)==0) next
        it <- iterm[1]
        inter_all[[length(inter_all)+1]] <- tibble(
          feature=ftr, covariate=cv,
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

    message("Done correlation: ", level_name, " | ", mname)
  }

  invisible(list(mat_f=mat_f, clr_f=clr_f, feat_anno=feat_anno))
}

## ---------- Run all levels ----------
# Host genus
lvl_host <- build_level("Host", mat_votu_in=mat_votu0, anno_votu_in=anno_votu0)
host_obj <- run_level_corr("Host", lvl_host)

# Host species
lvl_host_sp <- build_level("Host_species", mat_votu_in=mat_votu0, anno_votu_in=anno_votu0)
hostsp_obj <- run_level_corr("Host_species", lvl_host_sp)

# ViralTax ranks
for (rk in VIRAL_RANKS) {
  lvl_tax <- build_level("ViralTax", rank=rk, mat_votu_in=mat_votu0, anno_votu_in=anno_votu0)
  run_level_corr(paste0("ViralTax_", rk), lvl_tax)
}

## ---------- UPDRS scatterplots: Enterococcus phage ----------
if (MAKE_UPDRS_SCATTER) {

  # Prepare Host CLR matrices already computed in objects returned
  clr_host <- host_obj$clr_f           # features(host_genus) x samples
  clr_hostsp <- hostsp_obj$clr_f       # features(host_species_key) x samples
  feat_host <- host_obj$feat_anno
  feat_hostsp <- hostsp_obj$feat_anno

  scatter_root <- file.path(OUT_ROOT, "00_scatter_updrs3")
  dir_init(scatter_root)
  write_tsv0(meta, file.path(scatter_root, "00_manifest", "metadata_model_40.tsv"))

  ## Host_genus: Enterococcus
  if ("Enterococcus" %in% rownames(clr_host)) {
    for (m in MODELS) {
      out_dir <- file.path(scatter_root, "Host", m)
      dir_init(out_dir)

      x_raw <- as.numeric(clr_host["Enterococcus", sample_order])
      y_raw <- meta$updrs3_reduction_pct

      df_g <- make_xy_scatter(meta, x_raw, y_raw, model=m, mode="within_group")
      df_o <- make_xy_scatter(meta, x_raw, y_raw, model=m, mode="overall")

      write_tsv0(df_g |> mutate(level="Host_genus", feature="Enterococcus", model=m, plot_type="groupwise"),
                 file.path(out_dir, "tables", "scatter_data_groupwise.tsv"))
      write_tsv0(df_o |> mutate(level="Host_genus", feature="Enterococcus", model=m, plot_type="overall"),
                 file.path(out_dir, "tables", "scatter_data_overall.tsv"))

      xlab <- if (m=="unadjusted") "Enterococcus phage (CLR, Host_genus)" else paste0("Enterococcus phage (adj residual, ", m, ")")
      ylab <- if (m=="unadjusted") "MDS-UPDRS-III reduction (%)" else paste0("UPDRS-III reduction (adj residual, ", m, ")")

      plot_groupwise_scatter(df_g,
        title=paste0("Host_genus Enterococcus phage vs UPDRS-III reduction (", m, ")"),
        xlab=xlab, ylab=ylab,
        out_pdf=file.path(out_dir, "figures_pdf", "HostGenus_Enterococcus_groupwise.pdf"),
        out_png=file.path(out_dir, "figures_png", "HostGenus_Enterococcus_groupwise.png"))

      plot_overall_scatter(df_o,
        title=paste0("Host_genus Enterococcus phage vs UPDRS-III reduction (overall, ", m, ")"),
        xlab=xlab, ylab=ylab,
        out_pdf=file.path(out_dir, "figures_pdf", "HostGenus_Enterococcus_overall.pdf"),
        out_png=file.path(out_dir, "figures_png", "HostGenus_Enterococcus_overall.png"))
    }
  } else {
    writeLines("Host_genus feature 'Enterococcus' not found. Check host predictions / filtering.",
               con=file.path(scatter_root, "logs", "note_host_genus_missing.txt"))
  }

  ## Host_species: all Enterococcus species phage (features starting with 'Enterococcus ')
  enterococcus_sp_features <- rownames(clr_hostsp)[str_detect(rownames(clr_hostsp), "^Enterococcus\\b")]
  write_tsv0(tibble(host_species_feature = enterococcus_sp_features),
             file.path(scatter_root, "00_manifest", "enterococcus_host_species_features.tsv"))

  if (length(enterococcus_sp_features) > 0) {
    for (m in MODELS) {
      out_dir <- file.path(scatter_root, "Host_species", m)
      dir_init(out_dir)

      sum_rows <- list()
      long_g <- list()
      long_o <- list()

      for (sp in enterococcus_sp_features) {
        sp_slug <- safe_slug(sp)
        x_raw <- as.numeric(clr_hostsp[sp, sample_order])
        y_raw <- meta$updrs3_reduction_pct

        df_g <- make_xy_scatter(meta, x_raw, y_raw, model=m, mode="within_group")
        df_o <- make_xy_scatter(meta, x_raw, y_raw, model=m, mode="overall")

        long_g[[length(long_g)+1]] <- df_g |> mutate(level="Host_species", feature=sp, model=m, plot_type="groupwise")
        long_o[[length(long_o)+1]] <- df_o |> mutate(level="Host_species", feature=sp, model=m, plot_type="overall")

        st_gr <- spearman_fast(df_g$x_used[df_g$group=="PD_GR"], df_g$y_used[df_g$group=="PD_GR"])
        st_mr <- spearman_fast(df_g$x_used[df_g$group=="PD_MR"], df_g$y_used[df_g$group=="PD_MR"])
        st_all <- spearman_fast(df_o$x_used, df_o$y_used)

        sum_rows[[length(sum_rows)+1]] <- tibble(
          feature=sp, model=m,
          rho_GR=st_gr["rho"], p_GR=st_gr["p"], n_GR=st_gr["n"],
          rho_MR=st_mr["rho"], p_MR=st_mr["p"], n_MR=st_mr["n"],
          rho_all=st_all["rho"], p_all=st_all["p"], n_all=st_all["n"]
        )

        xlab <- if (m=="unadjusted") paste0(sp, " phage (CLR, Host_species)") else paste0(sp, " phage (adj residual, ", m, ")")
        ylab <- if (m=="unadjusted") "MDS-UPDRS-III reduction (%)" else paste0("UPDRS-III reduction (adj residual, ", m, ")")

        plot_groupwise_scatter(df_g,
          title=paste0("Host_species ", sp, " phage vs UPDRS-III reduction (", m, ")"),
          xlab=xlab, ylab=ylab,
          out_pdf=file.path(out_dir, "figures_pdf", paste0("HostSpecies_", sp_slug, "_groupwise.pdf")),
          out_png=file.path(out_dir, "figures_png", paste0("HostSpecies_", sp_slug, "_groupwise.png")))

        plot_overall_scatter(df_o,
          title=paste0("Host_species ", sp, " phage vs UPDRS-III reduction (overall, ", m, ")"),
          xlab=xlab, ylab=ylab,
          out_pdf=file.path(out_dir, "figures_pdf", paste0("HostSpecies_", sp_slug, "_overall.pdf")),
          out_png=file.path(out_dir, "figures_png", paste0("HostSpecies_", sp_slug, "_overall.png")))
      }

      write_tsv0(bind_rows(long_g), file.path(out_dir, "tables", "scatter_data_groupwise_all_long.tsv"))
      write_tsv0(bind_rows(long_o), file.path(out_dir, "tables", "scatter_data_overall_all_long.tsv"))
      write_tsv0(bind_rows(sum_rows), file.path(out_dir, "tables", "scatter_summary.tsv"))
    }
  } else {
    writeLines("No Host_species features starting with 'Enterococcus'. Check host predictions / parsing.",
               con=file.path(scatter_root, "logs", "note_host_species_missing.txt"))
  }
}

message("All virome correlations + Enterococcus UPDRS scatter done -> ", OUT_ROOT)