## =========================
## 12_multiomics_16S_virome_correlation_A_B_C_E_F.R
## Cross-omics correlation: 16S vs Virome
##
## Implements:
## (A) 16S Genus ↔ Virome Host_genus (overall only)
## (B) 16S Species ↔ Virome Host_species (overall only)
## (C) 16S rank ↔ Virome ViralTax rank (phylum/class/order/family/genus/species; overall only)
## (D) Heatmaps overall only (no group heatmaps)
## (E) Scatterplots (one-figure-per-file, your style):
##     - Genus: Enterococcus vs Host_genus Enterococcus phage
##     - Species: all Enterococcus spp. (16S) vs matching Host_species phage (intersection)
## (F) Bipartite network (bacteria—phage host) overall only, from (A)
##
## Models (3):
##   unadjusted: CLR
##   primary: residualize on group + age + bmi + disease_duration + constipation
##   sensitivity_add_dose: add levodopa_dose
##
## Inputs in wd:
##   metadata.csv
##   table_filtered_w_tax.xls   或   table_f_even_depth_w_tax.xls
##   vOTU_abundance_matrix_TPM.tsv
##   Virus_taxonomy.tsv (2 cols, no header)
##   cherry_prediction.tsv
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
INPUT_16S <- "table_f_even_depth_w_tax.xls"  #抽平/未抽平
TPM_FILE  <- "vOTU_abundance_matrix_TPM.tsv"
VTAX_FILE <- "Virus_taxonomy.tsv"
CHERRY_FILE <- "cherry_prediction.tsv"

## ---------- Global ----------
GROUP_LEVELS <- c("PD_GR","PD_MR")

MODELS <- c("unadjusted","primary","sensitivity_add_dose")

# 16S CLR
PSEUDOCOUNT_16S <- 0.5
# Virome CLR
PSEUDOCOUNT_VIROME <- 0.01

# Filtering (independent)
MIN_PREV_ALL <- 0.10      # >=4/40
MIN_TOTAL_16S <- 10
MIN_TOTAL_VIROME <- 0
MIN_PREV_BY_GROUP <- NA_real_  # optional: e.g. 0.05
MIN_MEAN_RA <- NA_real_        # optional: e.g. 1e-5

# Heatmap selection
P_CUTOFF <- 0.05
HEATMAP_MAX_ROWS <- 60
HEATMAP_MAX_COLS <- 60

# Network thresholds
NETWORK_Q_CUTOFF_MAIN <- 0.10
NETWORK_Q_CUTOFF_SENS <- 0.20

# Scatter colors (keep your style)
COL_MAP <- c(PD_GR = "#7CAE00", PD_MR = "#F8766D")
LAB_MAP <- c(PD_GR = "GR", PD_MR = "MR")
MAKE_PDF <- TRUE
MAKE_PNG <- TRUE

## ---------- Outputs ----------
OUT_ROOT <- "res_mo_even"  #抽平/未抽平
dir.create(OUT_ROOT, showWarnings = FALSE)

OUT_A <- file.path(OUT_ROOT, "A_Genus__HostGenus")
OUT_B <- file.path(OUT_ROOT, "B_Species__HostSpecies")
OUT_C <- file.path(OUT_ROOT, "C_LevelMatched")
OUT_E <- file.path(OUT_ROOT, "E_Scatter_Enterococcus")
OUT_F <- file.path(OUT_ROOT, "F_Network")

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

# Add suffix to virome feature names to distinguish from bacterial taxa
add_virus_suffix <- function(x) {
  x <- as.character(x)
  ifelse(stringr::str_ends(x, "_virus"), x, paste0(x, "_virus"))
}
strip_virus_suffix <- function(x) {
  sub("_virus$", "", as.character(x))
}

## ---------- Composition transforms ----------
clr_transform <- function(mat_feat_by_samp, pseudocount) {
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
  rowsum(mat, group=grp, reorder=TRUE)
}

## ---------- Residualization ----------
residualize_matrix <- function(Y_samp_by_feat, X_design) {
  qrX <- qr(X_design)
  qr.resid(qrX, Y_samp_by_feat)
}

build_design <- function(meta, model_name) {
  if (model_name == "unadjusted") return(NULL)
  if (model_name == "primary") {
    model.matrix(~ group + age + bmi + disease_duration + constipation, data = meta)
  } else {
    model.matrix(~ group + age + bmi + disease_duration + constipation + levodopa_dose, data = meta)
  }
}

## ---------- Spearman matrix (rank transform + Pearson) ----------
spearman_matrix <- function(X_feat_by_samp, Y_feat_by_samp) {
  stopifnot(ncol(X_feat_by_samp) == ncol(Y_feat_by_samp))

  n  <- ncol(X_feat_by_samp)
  nx <- nrow(X_feat_by_samp)
  ny <- nrow(Y_feat_by_samp)

  # 处理空矩阵/样本太少
  if (nx == 0 || ny == 0 || n < 6) {
    cor0 <- matrix(NA_real_, nrow = nx, ncol = ny,
                   dimnames = list(rownames(X_feat_by_samp), rownames(Y_feat_by_samp)))
    p0   <- cor0
    return(list(cor = cor0, p = p0, n = n))
  }

  rx <- t(apply(X_feat_by_samp, 1, rank, ties.method = "average"))
  ry <- t(apply(Y_feat_by_samp, 1, rank, ties.method = "average"))

  # cor() 在某些退化情形会返回向量/标量，这里强制成 nx×ny 矩阵
  cor_raw <- cor(t(rx), t(ry), method = "pearson")
  cor_mat <- matrix(as.numeric(cor_raw), nrow = nx, ncol = ny,
                    dimnames = list(rownames(X_feat_by_samp), rownames(Y_feat_by_samp)))

  r <- pmin(0.999999, pmax(-0.999999, cor_mat))
  tval <- r * sqrt((n - 2) / (1 - r^2))
  p_raw <- 2 * pt(-abs(tval), df = n - 2)
  p_mat <- matrix(as.numeric(p_raw), nrow = nx, ncol = ny, dimnames = dimnames(cor_mat))

  list(cor = cor_mat, p = p_mat, n = n)
}

## ---------- Heatmap plotting (overall only) ----------
plot_heatmap <- function(mat, out_pdf, out_png, title=NULL) {

  # >>> ADD: ensure output directories exist
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

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
    df <- as.data.frame(mat) |> rownames_to_column("row") |>
      pivot_longer(-row, names_to="col", values_to="val")
    p <- ggplot(df, aes(col, row, fill=val)) + geom_tile() + theme_bw() +
      theme(axis.text.x=element_text(angle=60, hjust=1)) + labs(x=NULL,y=NULL,title=title)
    if (MAKE_PDF) ggsave(out_pdf, p, width=10, height=max(4,0.12*nrow(mat)+2))
    if (MAKE_PNG) ggsave(out_png, p, width=10, height=max(4,0.12*nrow(mat)+2), dpi=300)
  }
}

## ---------- Scatter helpers (your style) ----------
spearman_fast <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 6) return(c(rho=NA_real_, p=NA_real_, n=n))
  rho <- suppressWarnings(cor(x, y, method="spearman"))
  rho2 <- min(0.999999, max(-0.999999, rho))
  tval <- rho2 * sqrt((n - 2) / (1 - rho2^2))
  p <- 2 * pt(-abs(tval), df = n - 2)
  c(rho=rho, p=p, n=n)
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

write_tsv0(meta, file.path(OUT_ROOT, "00_manifest_metadata_model_40.tsv"))

## ---------- Read 16S table and build rank matrices ----------
dt16 <- fread(INPUT_16S, sep="\t", header=TRUE, data.table=FALSE, check.names=FALSE)
id_col16 <- if ("#ID" %in% names(dt16)) "#ID" else names(dt16)[1]
tax_col16 <- if ("taxonomy" %in% names(dt16)) "taxonomy" else if ("taxon" %in% names(dt16)) "taxon" else NA_character_
stopifnot(!is.na(tax_col16))

sample_cols16 <- names(dt16)[str_detect(names(dt16), "^PD\\d+$")]
keep_samples16 <- intersect(sample_order, sample_cols16)
stopifnot(length(keep_samples16)==40)
keep_samples16 <- sample_order

asv_id <- as.character(dt16[[id_col16]])
mat_asv <- dt16[, keep_samples16, drop=FALSE]
mat_asv[] <- lapply(mat_asv, as.numeric)
mat_asv <- as.matrix(mat_asv)
rownames(mat_asv) <- asv_id

# parse 16S taxonomy ranks
parse_tax16 <- function(tax) {
  parts <- str_split(as.character(tax), ";", simplify=TRUE)
  get_rank <- function(row, key) {
    cand <- row[str_detect(row, paste0("^\\s*", key, "_+"))]
    if (length(cand)==0) return(NA_character_)
    v <- str_trim(cand[1])
    v <- str_replace(v, paste0("^\\s*", key, "_+"), "")
    if (v=="" || v=="NA") NA_character_ else v
  }
  tibble(
    phylum = apply(parts, 1, get_rank, key="p"),
    class  = apply(parts, 1, get_rank, key="c"),
    order  = apply(parts, 1, get_rank, key="o"),
    family = apply(parts, 1, get_rank, key="f"),
    genus  = apply(parts, 1, get_rank, key="g"),
    species= apply(parts, 1, get_rank, key="s")
  )
}
tax16_tbl <- tibble(asv_id=asv_id, taxonomy_raw=dt16[[tax_col16]]) |>
  bind_cols(parse_tax16(dt16[[tax_col16]])) |>
  mutate(across(c(phylum,class,order,family,genus,species), ~fill_unknown(.x, cur_column())))

write_tsv0(tax16_tbl, file.path(OUT_ROOT, "00_manifest_16S_taxonomy_parsed.tsv"))

# build 16S matrices for each rank
build_16S_rank <- function(rank_name) {
  map <- setNames(tax16_tbl[[rank_name]], tax16_tbl$asv_id)
  mat <- aggregate_sum(mat_asv, map)
  # filter
  f <- filter_matrix(mat, meta,
                     min_prev_all=MIN_PREV_ALL,
                     min_total=if (rank_name %in% c("genus","species")) MIN_TOTAL_16S else MIN_TOTAL_16S,
                     min_prev_by_group=MIN_PREV_BY_GROUP,
                     min_mean_ra=MIN_MEAN_RA)
  mat_f <- f$mat[, sample_order, drop=FALSE]
  # annotation: keep rank + parent ranks via mode
  ann <- tax16_tbl |>
    mutate(taxon = map[asv_id]) |>
    group_by(taxon) |>
    summarise(
      phylum = safe_mode(phylum),
      class  = safe_mode(class),
      order  = safe_mode(order),
      family = safe_mode(family),
      genus  = safe_mode(genus),
      species= safe_mode(species),
      n_asv = n(),
      .groups="drop"
    ) |>
    mutate(feature = taxon,
           level = paste0("16S_", rank_name)) |>
    filter(feature %in% rownames(mat_f))
  list(mat=mat_f, filter=f, anno=ann)
}

s16_phylum <- build_16S_rank("phylum")
s16_class  <- build_16S_rank("class")
s16_order  <- build_16S_rank("order")
s16_family <- build_16S_rank("family")
s16_genus  <- build_16S_rank("genus")
s16_species<- build_16S_rank("species")

## ---------- Read virome TPM and build Host/ViralTax matrices ----------
tpm_dt <- fread(TPM_FILE, sep="\t", header=TRUE, data.table=FALSE, check.names=FALSE)
stopifnot("id" %in% names(tpm_dt))
sample_cols_v <- setdiff(names(tpm_dt), "id")
sample_cols_v <- sample_cols_v[str_detect(sample_cols_v, "^PD\\d+$")]
keep_samples_v <- intersect(sample_order, sample_cols_v)
stopifnot(length(keep_samples_v)==40)
keep_samples_v <- sample_order

mat_votu <- tpm_dt[, keep_samples_v, drop=FALSE]
mat_votu[] <- lapply(mat_votu, as.numeric)
mat_votu <- as.matrix(mat_votu)
rownames(mat_votu) <- as.character(tpm_dt$id)

# viral taxonomy
vtax <- fread(VTAX_FILE, sep="\t", header=FALSE, data.table=FALSE)
colnames(vtax) <- c("votu_id","taxonomy_raw")
parse_vtax <- function(x) {
  parts <- str_split(as.character(x), ";", simplify=TRUE)
  get_rank <- function(row, prefix) {
    y <- row[str_detect(row, paste0("^\\s*", prefix))]
    if (length(y)==0) return(NA_character_)
    str_replace(str_trim(y[1]), paste0("^", prefix), "")
  }
  tibble(
    phylum  = apply(parts, 1, get_rank, prefix="p__"),
    class   = apply(parts, 1, get_rank, prefix="c__"),
    order   = apply(parts, 1, get_rank, prefix="o__"),
    family  = apply(parts, 1, get_rank, prefix="f__"),
    genus   = apply(parts, 1, get_rank, prefix="g__"),
    species = apply(parts, 1, get_rank, prefix="s__")
  )
}
vtax_tbl <- vtax |>
  mutate(votu_id=as.character(votu_id)) |>
  bind_cols(parse_vtax(vtax$taxonomy_raw)) |>
  mutate(across(c(phylum,class,order,family,genus,species), ~fill_unknown(.x, cur_column())))

# CHERRY host
ch <- fread(CHERRY_FILE, sep="\t", header=TRUE, data.table=FALSE)
# parse host genus/species
parse_host <- function(host_raw) {
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
  tibble(host_genus=host_genus, host_species=host_species)
}

host_tbl <- ch |>
  transmute(votu_id = as.character(Accession),
            Host_raw = as.character(Host),
            cherry_score = suppressWarnings(as.numeric(CHERRYScore)),
            method = as.character(Method))

host_tbl <- dplyr::bind_cols(host_tbl, parse_host(host_tbl$Host_raw)) |>
  mutate(host_genus = fill_unknown(host_genus, "host"),
         method = ifelse(is.na(method) | method=="", "Unspecified", method))

anno_votu <- tibble(votu_id=rownames(mat_votu)) |>
  left_join(vtax_tbl, by="votu_id") |>
  left_join(host_tbl, by="votu_id") |>
  mutate(host_genus = fill_unknown(host_genus, "host"))

# vOTU backbone filter (for consistent membership)
f0 <- filter_matrix(mat_votu, meta,
                    min_prev_all=MIN_PREV_ALL,
                    min_total=MIN_TOTAL_VIROME,
                    min_prev_by_group=MIN_PREV_BY_GROUP,
                    min_mean_ra=MIN_MEAN_RA)
mat_votu0 <- f0$mat[, sample_order, drop=FALSE]
anno_votu0 <- anno_votu |> filter(votu_id %in% rownames(mat_votu0))

# build virome Host_genus and Host_species matrices
build_virome_host_genus <- function() {
  map <- setNames(anno_votu0$host_genus, anno_votu0$votu_id)
  mat <- aggregate_sum(mat_votu0, map)
  f <- filter_matrix(mat, meta, MIN_PREV_ALL, MIN_TOTAL_VIROME, MIN_PREV_BY_GROUP, MIN_MEAN_RA)
  mat_f <- f$mat[, sample_order, drop=FALSE]

  # add suffix to virome feature names
  rownames(mat_f) <- add_virus_suffix(rownames(mat_f))
  f$mat <- mat_f
  f$diag <- f$diag |> dplyr::mutate(feature = add_virus_suffix(feature))

  ann <- anno_votu0 |>
    group_by(host_genus) |>
    summarise(n_votu=n(),
              member_votu_ids=paste(votu_id, collapse=";"),
              mean_score=if (all(is.na(cherry_score))) NA_real_ else mean(cherry_score, na.rm=TRUE),
              max_score=if (all(is.na(cherry_score))) NA_real_ else max(cherry_score, na.rm=TRUE),
              top_method=safe_mode(method),
              .groups="drop") |>
    mutate(
      feature_raw = host_genus,
      feature = add_virus_suffix(host_genus),
      level = "Virome_Host_genus"
    ) |>
    filter(feature %in% rownames(mat_f))

  list(mat=mat_f, filter=f, anno=ann)
}


build_virome_host_species <- function() {
  key <- ifelse(!is.na(anno_votu0$host_species) & anno_votu0$host_species!="",
                anno_votu0$host_species,
                paste0(anno_votu0$host_genus, " [genus]"))
  key <- fill_unknown(key, "host_species")
  map <- setNames(key, anno_votu0$votu_id)
  mat <- aggregate_sum(mat_votu0, map)
  f <- filter_matrix(mat, meta, MIN_PREV_ALL, MIN_TOTAL_VIROME, MIN_PREV_BY_GROUP, MIN_MEAN_RA)
  mat_f <- f$mat[, sample_order, drop=FALSE]

  # add suffix to virome feature names
  rownames(mat_f) <- add_virus_suffix(rownames(mat_f))
  f$mat <- mat_f
  f$diag <- f$diag |> dplyr::mutate(feature = add_virus_suffix(feature))

  ann <- tibble(votu_id=anno_votu0$votu_id,
                host_species_key=map[anno_votu0$votu_id],
                host_genus=anno_votu0$host_genus,
                cherry_score=anno_votu0$cherry_score,
                method=anno_votu0$method) |>
    group_by(host_species_key) |>
    summarise(n_votu=n(),
              member_votu_ids=paste(votu_id, collapse=";"),
              top_host_genus=safe_mode(host_genus),
              mean_score=if (all(is.na(cherry_score))) NA_real_ else mean(cherry_score, na.rm=TRUE),
              max_score=if (all(is.na(cherry_score))) NA_real_ else max(cherry_score, na.rm=TRUE),
              top_method=safe_mode(method),
              .groups="drop") |>
    mutate(
      feature_raw = host_species_key,
      feature = add_virus_suffix(host_species_key),
      level = "Virome_Host_species"
    ) |>
    filter(feature %in% rownames(mat_f))

  list(mat=mat_f, filter=f, anno=ann)
}


build_virome_vtax_rank <- function(rank_name) {
  map <- setNames(anno_votu0[[rank_name]], anno_votu0$votu_id)
  mat <- aggregate_sum(mat_votu0, map)
  f <- filter_matrix(mat, meta, MIN_PREV_ALL, MIN_TOTAL_VIROME, MIN_PREV_BY_GROUP, MIN_MEAN_RA)
  mat_f <- f$mat[, sample_order, drop=FALSE]

  # add suffix to virome feature names
  rownames(mat_f) <- add_virus_suffix(rownames(mat_f))
  f$mat <- mat_f
  f$diag <- f$diag |> dplyr::mutate(feature = add_virus_suffix(feature))

  ann <- anno_votu0 |>
    mutate(taxon = map[votu_id]) |>
    group_by(taxon) |>
    summarise(phylum=safe_mode(phylum),
              class=safe_mode(class),
              order=safe_mode(order),
              family=safe_mode(family),
              genus=safe_mode(genus),
              species=safe_mode(species),
              n_votu=n(),
              member_votu_ids=paste(votu_id, collapse=";"),
              top_host_genus=safe_mode(host_genus),
              top_host_species=safe_mode(host_species),
              .groups="drop") |>
    mutate(
      feature_raw = taxon,
      feature = add_virus_suffix(taxon),
      level = paste0("Virome_ViralTax_", rank_name)
    ) |>
    filter(feature %in% rownames(mat_f))

  list(mat=mat_f, filter=f, anno=ann)
}


v_host_genus <- build_virome_host_genus()
v_host_species <- build_virome_host_species()
v_vtax <- lapply(c("phylum","class","order","family","genus","species"), build_virome_vtax_rank)
names(v_vtax) <- c("phylum","class","order","family","genus","species")

## ---------- Core cross-omics runner ----------
run_pair_corr <- function(out_root, X_obj, Y_obj, pair_name) {
  dir_init(out_root)

  # save manifests
  write_tsv0(meta, file.path(out_root, "00_manifest", "metadata_model_40.tsv"))
  write_tsv0(X_obj$filter$summary, file.path(out_root, "00_manifest", "filter_summary_X.tsv"))
  write_tsv0(Y_obj$filter$summary, file.path(out_root, "00_manifest", "filter_summary_Y.tsv"))
  write_tsv0(X_obj$filter$diag, file.path(out_root, "00_manifest", "filter_diagnostics_X.tsv"))
  write_tsv0(Y_obj$filter$diag, file.path(out_root, "00_manifest", "filter_diagnostics_Y.tsv"))
  write_tsv0(X_obj$anno, file.path(out_root, "00_manifest", "X_feature_annotation.tsv"))
  write_tsv0(Y_obj$anno, file.path(out_root, "00_manifest", "Y_feature_annotation.tsv"))

  # compute CLR base
  X_clr <- clr_transform(X_obj$mat, pseudocount=PSEUDOCOUNT_16S)   # X features x samples
  Y_clr <- clr_transform(Y_obj$mat, pseudocount=PSEUDOCOUNT_VIROME)

  for (m in MODELS) {
    out_m <- file.path(out_root, m)
    dir_init(out_m)

    # complete-case samples for model covariates
    if (m == "unadjusted") {
      keep_idx <- rep(TRUE, nrow(meta))
    } else if (m == "primary") {
      keep_idx <- complete.cases(meta[, c("group","age","bmi","disease_duration","constipation"), drop=FALSE])
    } else {
      keep_idx <- complete.cases(meta[, c("group","age","bmi","disease_duration","constipation","levodopa_dose"), drop=FALSE])
    }
    keep_samples <- meta$sample_id[keep_idx]
    if (length(keep_samples) < 20) {
      writeLines(paste0("Too few complete-case samples for model ", m, ": n=", length(keep_samples)),
                 con=file.path(out_m, "logs", "note.txt"))
      next
    }

    X_use <- X_clr[, keep_samples, drop=FALSE]
    Y_use <- Y_clr[, keep_samples, drop=FALSE]

    # residualize for adjusted models
    if (m != "unadjusted") {
      Xdes <- build_design(meta[keep_idx, , drop=FALSE], m)
      # residualize samples x features
      Xres <- t(residualize_matrix(t(X_use), Xdes))
      Yres <- t(residualize_matrix(t(Y_use), Xdes))

      # 恢复行列名，防止后续写表/热图错位
      rownames(Xres) <- rownames(X_use); colnames(Xres) <- colnames(X_use)
      rownames(Yres) <- rownames(Y_use); colnames(Yres) <- colnames(Y_use)

      X_use <- Xres
      Y_use <- Yres
    }

    # save per-sample matrices used (reproducible)
    write_tsv0(as.data.frame(X_use) |> rownames_to_column("feature"),
               file.path(out_m, "per_sample", "X_values_used_features_by_sample.tsv"))
    write_tsv0(as.data.frame(Y_use) |> rownames_to_column("feature"),
               file.path(out_m, "per_sample", "Y_values_used_features_by_sample.tsv"))
    write_tsv0(meta[keep_idx, ], file.path(out_m, "per_sample", "sample_annotation_used.tsv"))

    # compute correlations
    cm <- spearman_matrix(X_use, Y_use)
    cor_mat <- cm$cor
    p_mat <- cm$p
    n_used <- cm$n

    # BH across all pairs
    p_vec <- as.vector(p_mat)
    q_vec <- p.adjust(p_vec, method="BH")
    q_mat <- matrix(q_vec, nrow=nrow(p_mat), ncol=ncol(p_mat),
                    dimnames=dimnames(p_mat))

    # long table with annotations
    long <- as.data.frame(as.table(cor_mat), stringsAsFactors=FALSE)
    colnames(long) <- c("X_feature","Y_feature","rho")
    long$p_value <- as.vector(p_mat)
    long$q_value <- as.vector(q_mat)
    long$n_used <- n_used

    long <- long |>
      left_join(X_obj$anno |> rename(X_feature=feature), by="X_feature") |>
      left_join(Y_obj$anno |> rename(Y_feature=feature), by="Y_feature")

    write_tsv0(long |> arrange(q_value),
               file.path(out_m, "tables", "pairwise_correlations_long.tsv"))

    # write full matrices (for reproducing heatmaps)
    write_tsv0(as.data.frame(cor_mat) |> rownames_to_column("X_feature"),
               file.path(out_m, "tables", "correlation_matrix_rho.tsv"))
    write_tsv0(as.data.frame(p_mat) |> rownames_to_column("X_feature"),
               file.path(out_m, "tables", "correlation_matrix_p.tsv"))
    write_tsv0(as.data.frame(q_mat) |> rownames_to_column("X_feature"),
               file.path(out_m, "tables", "correlation_matrix_q.tsv"))

    # heatmap subset selection (overall only)
    # pick top rows/cols by min p
    row_minp <- apply(p_mat, 1, function(x) suppressWarnings(min(x, na.rm=TRUE)))
    col_minp <- apply(p_mat, 2, function(x) suppressWarnings(min(x, na.rm=TRUE)))

    sel_rows <- names(sort(row_minp))[seq_len(min(HEATMAP_MAX_ROWS, length(row_minp)))]
    sel_cols <- names(sort(col_minp))[seq_len(min(HEATMAP_MAX_COLS, length(col_minp)))]

    # if there are p<0.05 rows/cols, prioritize them
    pr_rows <- names(row_minp)[is.finite(row_minp) & row_minp < P_CUTOFF]
    pr_cols <- names(col_minp)[is.finite(col_minp) & col_minp < P_CUTOFF]
    if (length(pr_rows) > 0) sel_rows <- unique(c(pr_rows, sel_rows))[seq_len(min(HEATMAP_MAX_ROWS, length(unique(c(pr_rows, sel_rows)))))]
    if (length(pr_cols) > 0) sel_cols <- unique(c(pr_cols, sel_cols))[seq_len(min(HEATMAP_MAX_COLS, length(unique(c(pr_cols, sel_cols)))))]
    sel_rows <- intersect(sel_rows, rownames(cor_mat))
    sel_cols <- intersect(sel_cols, colnames(cor_mat))

    sub_mat <- cor_mat[sel_rows, sel_cols, drop=FALSE]
    write_tsv0(as.data.frame(sub_mat) |> rownames_to_column("X_feature"),
               file.path(out_m, "tables", "heatmap_subset_rho.tsv"))

    plot_heatmap(sub_mat,
                 out_pdf=file.path(out_m, "figures_pdf", paste0("heatmap_rho_", pair_name, "_", m, ".pdf")),
                 out_png=file.path(out_m, "figures_png", paste0("heatmap_rho_", pair_name, "_", m, ".png")),
                 title=paste0(pair_name, " (", m, "), n=", n_used))

    writeLines(paste0("Done: ", pair_name, " | ", m, " | n_used=", n_used,
                      " | X=", nrow(X_use), " features, Y=", nrow(Y_use), " features"),
               con=file.path(out_m, "logs", "done.txt"))
  }
}

## ---------- Pack objects for A/B/C ----------
# 16S objects
X16_genus  <- list(mat=s16_genus$mat, filter=s16_genus$filter, anno=s16_genus$anno)
X16_species<- list(mat=s16_species$mat, filter=s16_species$filter, anno=s16_species$anno)
X16_rank <- list(phylum=s16_phylum, class=s16_class, order=s16_order, family=s16_family, genus=s16_genus, species=s16_species)

# virome objects
Yv_host_genus  <- list(mat=v_host_genus$mat, filter=v_host_genus$filter, anno=v_host_genus$anno)
Yv_host_species<- list(mat=v_host_species$mat, filter=v_host_species$filter, anno=v_host_species$anno)
Yv_rank <- v_vtax

## ---------- Run (A) ----------
run_pair_corr(OUT_A, X16_genus, Yv_host_genus, pair_name="16S_Genus__Virome_HostGenus")

## ---------- Run (B) ----------
run_pair_corr(OUT_B, X16_species, Yv_host_species, pair_name="16S_Species__Virome_HostSpecies")

## ---------- Run (C): level matched ----------
dir_init(OUT_C)
for (rk in c("phylum","class","order","family","genus","species")) {
  out_rk <- file.path(OUT_C, paste0("rank_", rk))
  run_pair_corr(out_rk,
                list(mat=X16_rank[[rk]]$mat, filter=X16_rank[[rk]]$filter, anno=X16_rank[[rk]]$anno),
                list(mat=Yv_rank[[rk]]$mat, filter=Yv_rank[[rk]]$filter, anno=Yv_rank[[rk]]$anno),
                pair_name=paste0("16S_", rk, "__ViromeViralTax_", rk))
}

## ---------- (E) Scatterplots ----------
dir_init(OUT_E)

# Helper to get adjusted x/y for scatter under each model
get_xy_by_model <- function(x_vec, y_vec, model_name, meta_df) {
  if (model_name == "unadjusted") return(list(x=x_vec, y=y_vec, keep=rep(TRUE, nrow(meta_df))))
  if (model_name == "primary") {
    keep <- complete.cases(meta_df[, c("group","age","bmi","disease_duration","constipation"), drop=FALSE])
    Xdes <- build_design(meta_df[keep, , drop=FALSE], "primary")
  } else {
    keep <- complete.cases(meta_df[, c("group","age","bmi","disease_duration","constipation","levodopa_dose"), drop=FALSE])
    Xdes <- build_design(meta_df[keep, , drop=FALSE], "sensitivity_add_dose")
  }
  # residualize x and y using same design matrix
  xr <- rep(NA_real_, length(x_vec)); yr <- rep(NA_real_, length(y_vec))
  fitx <- lm.fit(Xdes, x_vec[keep]); fity <- lm.fit(Xdes, y_vec[keep])
  xr[keep] <- fitx$residuals
  yr[keep] <- fity$residuals
  list(x=xr, y=yr, keep=keep)
}

# (E1) Genus: Enterococcus vs Enterococcus phage host_genus
if ("Enterococcus" %in% rownames(s16_genus$mat) && add_virus_suffix("Enterococcus") %in% rownames(v_host_genus$mat)) {
  for (m in MODELS) {
    out_m <- file.path(OUT_E, "Genus_Enterococcus__HostGenus_Enterococcus", m)
    dir_init(out_m)

   x_raw <- as.numeric(clr_transform(s16_genus$mat, PSEUDOCOUNT_16S)["Enterococcus", sample_order])
    y_raw <- as.numeric(clr_transform(v_host_genus$mat, PSEUDOCOUNT_VIROME)[add_virus_suffix("Enterococcus"), sample_order])

    xy <- get_xy_by_model(x_raw, y_raw, m, meta)
    df <- meta |> mutate(x_used=xy$x, y_used=xy$y)

    write_tsv0(df |> mutate(model=m),
               file.path(out_m, "tables", "scatter_data.tsv"))

    xlab <- if (m=="unadjusted") "16S Enterococcus (CLR)" else paste0("16S Enterococcus (resid, ", m, ")")
    ylab <- if (m=="unadjusted") "Enterococcus phage (Host_genus CLR)" else paste0("Enterococcus phage (resid, ", m, ")")

    plot_groupwise_scatter(df,
      title=paste0("16S Enterococcus vs Enterococcus phage (groupwise, ", m, ")"),
      xlab=xlab, ylab=ylab,
      out_pdf=file.path(out_m, "figures_pdf", "scatter_groupwise.pdf"),
      out_png=file.path(out_m, "figures_png", "scatter_groupwise.png"))

    plot_overall_scatter(df,
      title=paste0("16S Enterococcus vs Enterococcus phage (overall, ", m, ")"),
      xlab=xlab, ylab=ylab,
      out_pdf=file.path(out_m, "figures_pdf", "scatter_overall.pdf"),
      out_png=file.path(out_m, "figures_png", "scatter_overall.png"))
  }
} else {
  writeLines("Enterococcus not found in either 16S genus or virome host_genus matrices after filtering.",
             con=file.path(OUT_E, "logs", "note_enterococcus_missing_genus.txt"))
}

# (E2) Species: all Enterococcus spp intersection
s16_sp_names <- rownames(s16_species$mat)[str_detect(rownames(s16_species$mat), "^Enterococcus\\b")]
v_sp_names <- strip_virus_suffix(rownames(v_host_species$mat)[str_detect(strip_virus_suffix(rownames(v_host_species$mat)), "^Enterococcus\b")])
sp_intersect <- intersect(s16_sp_names, v_sp_names)

cand_long <- dplyr::bind_rows(
  tibble(source = "16S_species", species = s16_sp_names),
  tibble(source = "Host_species", species = v_sp_names)
)

write_tsv0(cand_long,
           file.path(OUT_E, "00_manifest", "Enterococcus_species_candidates.tsv"))
write_tsv0(tibble(species_intersection = sp_intersect),
           file.path(OUT_E, "00_manifest", "Enterococcus_species_intersection.tsv"))

if (length(sp_intersect) > 0) {
  for (m in MODELS) {
    out_m <- file.path(OUT_E, "Species_Enterococcus_pairs", m)
    dir_init(out_m)

    # precompute CLR matrices for speed
    clr16_sp <- clr_transform(s16_species$mat, PSEUDOCOUNT_16S)
    clrv_sp  <- clr_transform(v_host_species$mat, PSEUDOCOUNT_VIROME)

    sum_rows <- list()

    for (sp in sp_intersect) {
      sp_slug <- safe_slug(sp)

      x_raw <- as.numeric(clr16_sp[sp, sample_order])
      y_raw <- as.numeric(clrv_sp[add_virus_suffix(sp), sample_order])

      xy <- get_xy_by_model(x_raw, y_raw, m, meta)
      df <- meta |> mutate(x_used=xy$x, y_used=xy$y)

      write_tsv0(df |> mutate(model=m, species=sp),
                 file.path(out_m, "tables", paste0("scatter_data_", sp_slug, ".tsv")))

      xlab <- if (m=="unadjusted") paste0("16S ", sp, " (CLR)") else paste0("16S ", sp, " (resid, ", m, ")")
      ylab <- if (m=="unadjusted") paste0("Phage host ", sp, " (CLR)") else paste0("Phage host ", sp, " (resid, ", m, ")")

      plot_groupwise_scatter(df,
        title=paste0("16S ", sp, " vs phage host ", sp, " (groupwise, ", m, ")"),
        xlab=xlab, ylab=ylab,
        out_pdf=file.path(out_m, "figures_pdf", paste0("scatter_groupwise_", sp_slug, ".pdf")),
        out_png=file.path(out_m, "figures_png", paste0("scatter_groupwise_", sp_slug, ".png")))

      plot_overall_scatter(df,
        title=paste0("16S ", sp, " vs phage host ", sp, " (overall, ", m, ")"),
        xlab=xlab, ylab=ylab,
        out_pdf=file.path(out_m, "figures_pdf", paste0("scatter_overall_", sp_slug, ".pdf")),
        out_png=file.path(out_m, "figures_png", paste0("scatter_overall_", sp_slug, ".png")))

      st_all <- spearman_fast(df$x_used, df$y_used)
      sum_rows[[length(sum_rows)+1]] <- tibble(species=sp, model=m, rho_all=st_all["rho"], p_all=st_all["p"], n_all=st_all["n"])
    }

    write_tsv0(bind_rows(sum_rows) |> arrange(p_all),
               file.path(out_m, "tables", "summary_all_species_overall.tsv"))
  }
} else {
  writeLines("No intersection species between 16S Enterococcus species and virome Host_species Enterococcus species.",
             con=file.path(OUT_E, "logs", "note_enterococcus_species_no_intersection.txt"))
}

## ---------- (F) Bipartite network from (A) ----------
dir_init(OUT_F)

make_network <- function(edge_df, out_prefix) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    writeLines("igraph not installed; network plot skipped. Edge list still written.",
               con = paste0(out_prefix, "_note.txt"))
    return(invisible(NULL))
  }

  nodes_b <- unique(edge_df$X_feature)
  nodes_v <- unique(edge_df$Y_feature)
  nodes <- tibble(
    name = c(nodes_b, nodes_v),
    type = c(rep("bacteria", length(nodes_b)), rep("phage_host", length(nodes_v)))
  )

  edges <- edge_df |>
    transmute(
      from = X_feature,
      to   = Y_feature,
      rho  = as.numeric(rho),
      weight_pos = pmax(1e-6, abs(as.numeric(rho)))
    ) |>
    filter(is.finite(rho))

  if (nrow(edges) < 2) {
    writeLines("Too few edges to plot network.",
               con = paste0(out_prefix, "_note_too_few_edges.txt"))
    return(invisible(NULL))
  }

  g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  igraph::E(g)$rho <- as.numeric(igraph::E(g)$rho)
  igraph::E(g)$weight_pos <- as.numeric(igraph::E(g)$weight_pos)

  set.seed(123)
  lay <- igraph::layout_with_fr(g, weights = igraph::E(g)$weight_pos)

  plot_one <- function() {
    plot(g,
         vertex.label.cex = 0.6,
         vertex.size = 5,
         edge.width = pmax(1, igraph::E(g)$weight_pos * 3),
         layout = lay)
  }

  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

  if (MAKE_PDF) {
    pdf(paste0(out_prefix, ".pdf"), width = 10, height = 7)
    plot_one()
    dev.off()
  }
  if (MAKE_PNG) {
    png(paste0(out_prefix, ".png"), width = 2600, height = 1800, res = 300)
    plot_one()
    dev.off()
  }

  invisible(NULL)
}

# Network based on A results (each model)
for (m in MODELS) {
  fpath <- file.path(OUT_A, m, "tables", "pairwise_correlations_long.tsv")
  if (!file.exists(fpath)) next
  ed <- readr::read_tsv(fpath, show_col_types=FALSE)

  # main threshold
  ed_main <- ed |> filter(!is.na(q_value) & q_value < NETWORK_Q_CUTOFF_MAIN) |> arrange(q_value)
  write_tsv0(ed, file.path(OUT_F, paste0("edges_full_", m, ".tsv")))
  write_tsv0(ed_main, file.path(OUT_F, paste0("edges_q_lt_", NETWORK_Q_CUTOFF_MAIN, "_", m, ".tsv")))

  # sensitivity threshold
  ed_sens <- ed |> filter(!is.na(q_value) & q_value < NETWORK_Q_CUTOFF_SENS) |> arrange(q_value)
  write_tsv0(ed_sens, file.path(OUT_F, paste0("edges_q_lt_", NETWORK_Q_CUTOFF_SENS, "_", m, ".tsv")))

  # plot
  if (nrow(ed_main) >= 2) {
    make_network(ed_main, file.path(OUT_F, paste0("network_q", NETWORK_Q_CUTOFF_MAIN, "_", m)))
  } else if (nrow(ed_sens) >= 2) {
    make_network(ed_sens, file.path(OUT_F, paste0("network_q", NETWORK_Q_CUTOFF_SENS, "_", m)))
  } else {
    writeLines(paste0("Too few edges for network under q<", NETWORK_Q_CUTOFF_MAIN, " or q<", NETWORK_Q_CUTOFF_SENS, " in model ", m),
               con=file.path(OUT_F, paste0("note_too_few_edges_", m, ".txt")))
  }
}

message("All multi-omics correlation analyses done -> ", OUT_ROOT)