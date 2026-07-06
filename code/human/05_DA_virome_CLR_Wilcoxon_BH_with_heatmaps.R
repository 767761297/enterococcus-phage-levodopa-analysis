## =========================
## Virome DA + Heatmaps (TaxHost full, 3 models, CLR+Wilcoxon+BH)
## - Levels: Host(genus) + Host_species + ViralTax_(phylum/class/order/family/genus/species) + vOTU(stats-only)
## - Heatmaps: Host + Host_species + each ViralTax rank (NO vOTU heatmap)
## - Tables: full stats + full annotations + per-sample ALL (wide as requested) + selected
## - All figures output PDF (PNG optional)
## Inputs (in working dir):
##   metadata.csv
##   vOTU_abundance_matrix_TPM.tsv
##   Virus_taxonomy.tsv          (2 cols, no header)
##   cherry_prediction.tsv
## Output root:
##   results_virome_DA_CLR_Wilcoxon_BH/
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

## ---------- inputs ----------
META_FILE <- "metadata.csv"
TPM_FILE  <- "vOTU_abundance_matrix_TPM.tsv"
TAX_FILE  <- "Virus_taxonomy.tsv"       # 2 cols, no header
HOST_FILE <- "cherry_prediction.tsv"

## ---------- output ----------
OUT_ROOT <- "results_virome_DA_CLR_Wilcoxon_BH"
dir.create(OUT_ROOT, showWarnings = FALSE)

## ---------- global ----------
GROUP_LEVELS <- c("PD_GR","PD_MR")   # GR first; effect = MR - GR
FDR_CUTOFF <- 0.10

## Filtering (independent; does not use group labels)
MIN_PREVALENCE_ALL <- 0.10
MIN_TOTAL_TPM <- 0

# optional stricter independent filters (set NA to disable)
MIN_PREV_BY_GROUP <- NA_real_   # e.g. 0.05
MIN_MEAN_RA <- NA_real_         # e.g. 1e-5

## CLR
PSEUDOCOUNT <- 0.01

## Heatmap inclusion (per level/model): include p<0.05; if none, fallback top N by p
P_CUTOFF <- 0.05
FALLBACK_TOP_N <- 30

## Heatmap plot settings
CLUSTER_ROWS_DEFAULT <- TRUE
CLUSTER_COLS_DEFAULT <- FALSE
MAKE_PDF <- TRUE
MAKE_PNG <- TRUE

## Which taxonomy ranks to run
VIRAL_RANKS <- c("phylum","class","order","family","genus","species")

## Heatmaps requested
MAKE_HEATMAP_HOST <- TRUE
MAKE_HEATMAP_TAX  <- TRUE

## vOTU-level: stats + selected tables only (no heatmap)
DO_vOTU_LEVEL <- TRUE
vOTU_EXPORT_PER_SAMPLE_ALL <- FALSE   # huge if TRUE; default FALSE

## Presence supplement (vOTU-level)
RUN_vOTU_PRESENCE_SUPPLEMENT <- TRUE
PRES_MIN_PREVALENCE_ALL <- 0.10
PRES_TEST_METHOD <- "auto"  # auto/fisher/chisq

## Membership ID collapsing (to avoid gigantic cells)
MAX_MEMBER_IDS_IN_CELL <- 200

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
write_feature_sample_wide_with_stats <- function(test_feat_by_samp,
                                                 res_df,
                                                 sample_order,
                                                 out_all_path,
                                                 out_p05_path,
                                                 p_cutoff = 0.05) {
  # test_feat_by_samp: features x samples (colnames are sample IDs)
  # res_df must contain: feature, p_value, q_value, effect_median_diff, rank_biserial, W, U, direction

  # 1) feature x sample (wide)
  wide_vals <- as.data.frame(test_feat_by_samp[, sample_order, drop = FALSE]) |>
    tibble::rownames_to_column("feature")

  # 2) stats columns
  stats_cols <- res_df |>
    select(feature,
           effect_median_diff, rank_biserial, W, U,
           p_value, q_value, direction)

  # 3) merge and order columns: feature + sample cols + stats cols
  wide_all <- wide_vals |>
    left_join(stats_cols, by = "feature")

  wide_all <- wide_all[, c("feature", sample_order,
                           "effect_median_diff", "rank_biserial", "W", "U",
                           "p_value", "q_value", "direction"), drop = FALSE]

  write_tsv0(wide_all, out_all_path)

  # 4) p<0.05 filtered
  wide_p <- wide_all |>
    filter(!is.na(p_value) & p_value < p_cutoff)

  write_tsv0(wide_p, out_p05_path)
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

collapse_ids <- function(ids, max_n = 200) {
  ids <- as.character(ids)
  ids <- ids[!is.na(ids) & ids != ""]
  if (length(ids) == 0) return(NA_character_)
  if (length(ids) <= max_n) return(paste(ids, collapse = ";"))
  head_part <- paste(ids[seq_len(max_n)], collapse = ";")
  paste0(head_part, ";...(+", length(ids) - max_n, " more)")
}

row_zscore <- function(mat_feat_by_samp) {
  mu <- rowMeans(mat_feat_by_samp, na.rm = TRUE)
  sdv <- apply(mat_feat_by_samp, 1, sd, na.rm = TRUE)
  sdv[sdv == 0 | is.na(sdv)] <- 1
  sweep(sweep(mat_feat_by_samp, 1, mu, "-"), 1, sdv, "/")
}

clr_transform <- function(mat_feat_by_samp, pseudocount = 0.01) {
  x <- mat_feat_by_samp + pseudocount
  x <- sweep(x, 2, colSums(x), "/")
  lx <- log(x)
  sweep(lx, 2, colMeans(lx), "-")
}

residualize_matrix <- function(Y_samp_by_feat, X_design) {
  qrX <- qr(X_design)
  qr.resid(qrX, Y_samp_by_feat)
}

# Wilcoxon summary: U from ranks => rank-biserial in [-1,1]
wilcox_one <- function(x_mr, x_gr) {
  n1 <- length(x_mr); n2 <- length(x_gr)
  wt <- suppressWarnings(wilcox.test(x_mr, x_gr, exact = FALSE))
  p <- wt$p.value
  Wstat <- as.numeric(wt$statistic)

  z <- c(x_mr, x_gr)
  r <- rank(z, ties.method = "average")
  Wrank <- sum(r[seq_len(n1)])
  U <- Wrank - n1*(n1+1)/2
  rb <- (2*U)/(n1*n2) - 1
  rb <- max(-1, min(1, rb))
  med_diff <- median(x_mr, na.rm = TRUE) - median(x_gr, na.rm = TRUE)
  c(p_value=p, W=Wstat, U=U, rank_biserial=rb, median_diff=med_diff)
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
    diagnostics = tibble(feature=rownames(mat_feat_by_samp),
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

# NOW outputs host_genus + host_species; host_display stays genus (for Host-level)
parse_host_display <- function(host_field) {
  host_field <- as.character(host_field)
  host_field[host_field %in% c("-", "", "NA")] <- NA_character_

  rank <- ifelse(!is.na(host_field) & str_detect(host_field, ":"), str_split_fixed(host_field, ":", 2)[,1], NA_character_)
  name <- ifelse(!is.na(host_field) & str_detect(host_field, ":"), str_split_fixed(host_field, ":", 2)[,2], host_field)
  name <- str_trim(name)

  host_species <- ifelse(!is.na(rank) & rank=="species" & !is.na(name) & name!="", name, NA_character_)
  host_genus <- ifelse(!is.na(rank) & rank=="genus" & !is.na(name) & name!="", name,
                       ifelse(!is.na(host_species) & str_detect(host_species, "\\s+"),
                              str_split_fixed(host_species, "\\s+", 2)[,1],
                              NA_character_))
  host_genus <- ifelse(!is.na(host_genus) & host_genus!="", host_genus, "Unassigned_host")

  tibble(
    host_rank   = rank,
    host_name   = name,
    host_genus  = host_genus,
    host_species= host_species,
    host_display= host_genus
  )
}

plot_heatmap <- function(mat_z, anno_col, out_pdf, out_png,
                         cluster_rows = TRUE, cluster_cols = FALSE) {
  # disable clustering if dimension <2
  cr <- isTRUE(cluster_rows) && nrow(mat_z) >= 2
  cc <- isTRUE(cluster_cols) && ncol(mat_z) >= 2

  if (requireNamespace("pheatmap", quietly = TRUE)) {
    ann <- as.data.frame(anno_col); rownames(ann) <- colnames(mat_z)
    if (MAKE_PDF) {
      pdf(out_pdf, width = 9, height = max(4, 0.18*nrow(mat_z)+2))
      pheatmap::pheatmap(mat_z,
                         cluster_rows = cr,
                         cluster_cols = cc,
                         annotation_col = ann,
                         show_colnames = FALSE,
                         fontsize_row = 8)
      dev.off()
    }
    if (MAKE_PNG) {
      png(out_png, width = 2200, height = max(1400, 40*nrow(mat_z)+600), res = 300)
      pheatmap::pheatmap(mat_z,
                         cluster_rows = cr,
                         cluster_cols = cc,
                         annotation_col = ann,
                         show_colnames = FALSE,
                         fontsize_row = 8)
      dev.off()
    }
  } else {
    df <- as.data.frame(mat_z) |>
      rownames_to_column("feature") |>
      pivot_longer(-feature, names_to="sample_id", values_to="z")
    df$sample_id <- factor(df$sample_id, levels = colnames(mat_z))
    df$feature <- factor(df$feature, levels = rev(rownames(mat_z)))

    p <- ggplot(df, aes(sample_id, feature, fill=z)) +
      geom_tile() +
      theme_bw() +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      labs(x=NULL, y=NULL, fill="Z")

    if (MAKE_PDF) ggsave(out_pdf, p, width=9, height=max(4,0.18*nrow(mat_z)+2))
    if (MAKE_PNG) ggsave(out_png, p, width=9, height=max(4,0.18*nrow(mat_z)+2), dpi=300)
  }
}

## wide writer: rows = samples, cols = features
write_matrix_wide_sample_rows <- function(mat_feat_by_samp, out_path, sample_order) {
  # mat_feat_by_samp: features x samples
  mat <- mat_feat_by_samp[, sample_order, drop=FALSE]
  wide <- as.data.frame(t(mat))
  wide <- tibble::rownames_to_column(wide, "sample_id")
  write_tsv0(wide, out_path)
}

## ---------- load metadata ----------
meta_raw <- read_csv_smart(META_FILE) |> janitor::clean_names(replace=c("\u03bc"="u"))
sample_col <- find_col(meta_raw, c("^sample$"))
group_col  <- find_col(meta_raw, c("^group$"))
age_col    <- find_col(meta_raw, c("^age$"))
bmi_col    <- find_col(meta_raw, c("^bmi"))
dur_col    <- find_col(meta_raw, c("disease_duration"))
con_col    <- find_col(meta_raw, c("constipation"))
dose_col   <- find_col(meta_raw, c("levodopa.*dose|levodopa_daily_dose"))

stopifnot(!is.na(sample_col), !is.na(group_col))

meta <- meta_raw |>
  mutate(sample_id = normalize_id(.data[[sample_col]]),
         group = .data[[group_col]]) |>
  filter(!is.na(sample_id) & sample_id!="") |>
  filter(group %in% GROUP_LEVELS) |>
  transmute(
    sample_id,
    group = factor(group, levels = GROUP_LEVELS),
    age = as.numeric(.data[[age_col]]),
    bmi = as.numeric(.data[[bmi_col]]),
    disease_duration = as.numeric(.data[[dur_col]]),
    constipation = factor(as.integer(.data[[con_col]]), levels=c(0,1), labels=c("No","Yes")),
    levodopa_dose = if(!is.na(dose_col)) suppressWarnings(as.numeric(.data[[dose_col]])) else NA_real_
  )

stopifnot(nrow(meta)==40)

# order samples: GR then MR, within group by sample_id
meta <- meta |> arrange(group, sample_id)
sample_order <- meta$sample_id

## ---------- read TPM ----------
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

## ---------- read annotations ----------
tax_dt <- fread(TAX_FILE, sep="\t", header=FALSE, data.table=FALSE)
colnames(tax_dt) <- c("votu_id","taxonomy_raw")
tax_tbl <- tax_dt |>
  mutate(votu_id = as.character(votu_id)) |>
  bind_cols(parse_virus_tax(tax_dt$taxonomy_raw)) |>
  mutate(across(c(realm,kingdom,phylum,class,order,family,genus,species),
                ~fill_unknown(.x, cur_column())))

host_dt <- fread(HOST_FILE, sep="\t", header=TRUE, data.table=FALSE)
host_tbl <- host_dt |>
  transmute(
    votu_id = as.character(Accession),
    length = suppressWarnings(as.numeric(Length)),
    Host_raw = as.character(Host),
    cherry_score = suppressWarnings(as.numeric(CHERRYScore)),
    method = as.character(Method),
    host_ncbi_lineage_raw = as.character(Host_NCBI_lineage),
    host_gtdb_lineage_raw = as.character(Host_GTDB_lineage)
  )
host_tbl <- bind_cols(host_tbl, parse_host_display(host_tbl$Host_raw)) |>
  mutate(
    host_genus = fill_unknown(host_genus, "host"),
    host_display = host_genus,
    method = ifelse(is.na(method) | method=="", "Unspecified", method)
  )

anno_votu <- tibble(votu_id = rownames(mat_votu)) |>
  left_join(tax_tbl, by="votu_id") |>
  left_join(host_tbl, by="votu_id") |>
  mutate(
    host_genus = fill_unknown(host_genus, "host"),
    host_display = host_genus
  )

## ---------- manifest ----------
dir_init(file.path(OUT_ROOT, "00_manifest"))
write_tsv0(meta, file.path(OUT_ROOT, "00_manifest", "metadata_model_40.tsv"))
write_tsv0(anno_votu, file.path(OUT_ROOT, "00_manifest", "annotation_vOTU_full.tsv"))
write_tsv0(tibble(
  MIN_PREVALENCE_ALL=MIN_PREVALENCE_ALL,
  MIN_TOTAL_TPM=MIN_TOTAL_TPM,
  MIN_PREV_BY_GROUP=MIN_PREV_BY_GROUP,
  MIN_MEAN_RA=MIN_MEAN_RA,
  PSEUDOCOUNT=PSEUDOCOUNT,
  P_CUTOFF=P_CUTOFF,
  FALLBACK_TOP_N=FALLBACK_TOP_N,
  FDR_CUTOFF=FDR_CUTOFF,
  MAX_MEMBER_IDS_IN_CELL=MAX_MEMBER_IDS_IN_CELL
), file.path(OUT_ROOT, "00_manifest", "params.tsv"))

## ---------- models ----------
X_primary <- model.matrix(~ age + bmi + disease_duration + constipation, data = meta)
X_sens    <- model.matrix(~ age + bmi + disease_duration + constipation + levodopa_dose, data = meta)

MODELS <- list(
  unadjusted = list(type="clr", X=NULL),
  primary = list(type="resid", X=X_primary),
  sensitivity_add_dose = list(type="resid", X=X_sens)
)

## ---------- base vOTU filtering (common backbone across levels) ----------
f0 <- filter_matrix(mat_votu, meta,
                    min_prev_all = MIN_PREVALENCE_ALL,
                    min_total = MIN_TOTAL_TPM,
                    min_prev_by_group = MIN_PREV_BY_GROUP,
                    min_mean_ra = MIN_MEAN_RA)
mat_votu0 <- f0$mat
anno_votu0 <- anno_votu |> filter(votu_id %in% rownames(mat_votu0))
write_tsv0(f0$summary, file.path(OUT_ROOT, "00_manifest", "filter_summary_vOTU_base.tsv"))
saveRDS(mat_votu0, file.path(OUT_ROOT, "00_manifest", "TPM_vOTU_base_filtered.rds"), compress="xz")
write_tsv0(anno_votu0, file.path(OUT_ROOT, "00_manifest", "annotation_vOTU_base_filtered.tsv"))

## ---------- build levels ----------
build_level <- function(level, rank = NULL, mat_votu_in, anno_votu_in) {

  if (level == "vOTU") {
    mat <- mat_votu_in
    feat_anno <- anno_votu_in |>
      rename(feature = votu_id) |>
      mutate(
        level = "vOTU",
        display_label = paste0(host_display, " [", feature, "]")
      )
    map_tbl <- tibble(feature = feat_anno$feature, member_votu_id = feat_anno$feature)
    return(list(mat=mat, feat_anno=feat_anno, map_tbl=map_tbl))
  }

  if (level == "Host") {
    map <- setNames(anno_votu_in$host_genus, anno_votu_in$votu_id)
    mat <- aggregate_sum(mat_votu_in, map)

    host_sum <- anno_votu_in |>
      group_by(host_genus) |>
      summarise(
        n_member_votu = n(),
        member_votu_ids = collapse_ids(votu_id, MAX_MEMBER_IDS_IN_CELL),
        # viral taxonomy summaries
        top_phylum = safe_mode(phylum),
        top_class  = safe_mode(class),
        top_order  = safe_mode(order),
        top_family = safe_mode(family),
        # species info retained as summary
        n_distinct_species = dplyr::n_distinct(host_species[!is.na(host_species)]),
        top_host_species = safe_mode(host_species),
        # scores
        mean_cherry_score = if (all(is.na(cherry_score))) NA_real_ else mean(cherry_score, na.rm=TRUE),
        max_cherry_score  = if (all(is.na(cherry_score))) NA_real_ else max(cherry_score, na.rm=TRUE),
        top_method = safe_mode(method),
        .groups="drop"
      ) |>
      mutate(
        feature = host_genus,
        level = "Host_genus",
        display_label = ifelse(!is.na(top_host_species) & top_host_species!="",
                               paste0(host_genus, " | top:", top_host_species),
                               host_genus)
      ) |>
      select(feature, display_label, level,
             n_member_votu, member_votu_ids,
             top_phylum, top_class, top_order, top_family,
             n_distinct_species, top_host_species,
             mean_cherry_score, max_cherry_score, top_method)

    map_tbl <- tibble(votu_id = names(map), host_genus = unname(map))
    return(list(mat=mat, feat_anno=host_sum, map_tbl=map_tbl))
  }

  if (level == "Host_species") {
    host_species_key <- ifelse(!is.na(anno_votu_in$host_species) & anno_votu_in$host_species != "",
                               anno_votu_in$host_species,
                               paste0(anno_votu_in$host_genus, " [genus]"))
    host_species_key <- fill_unknown(host_species_key, "host_species")

    map <- setNames(host_species_key, anno_votu_in$votu_id)
    mat <- aggregate_sum(mat_votu_in, map)

    sp_sum <- anno_votu_in |>
      mutate(host_species_key = map[votu_id]) |>
      group_by(host_species_key) |>
      summarise(
        n_member_votu = n(),
        member_votu_ids = collapse_ids(votu_id, MAX_MEMBER_IDS_IN_CELL),
        top_host_genus = safe_mode(host_genus),
        mean_cherry_score = if (all(is.na(cherry_score))) NA_real_ else mean(cherry_score, na.rm=TRUE),
        max_cherry_score  = if (all(is.na(cherry_score))) NA_real_ else max(cherry_score, na.rm=TRUE),
        top_family = safe_mode(family),
        top_method = safe_mode(method),
        .groups="drop"
      ) |>
      mutate(
        feature = host_species_key,
        level = "Host_species",
        display_label = host_species_key
      ) |>
      select(feature, display_label, level,
             n_member_votu, member_votu_ids,
             top_host_genus, top_family,
             mean_cherry_score, max_cherry_score, top_method)

    map_tbl <- tibble(votu_id = names(map), host_species = unname(map))
    return(list(mat=mat, feat_anno=sp_sum, map_tbl=map_tbl))
  }

  if (level == "ViralTax") {
    stopifnot(!is.null(rank))
    taxon <- anno_votu_in[[rank]]
    taxon <- fill_unknown(taxon, rank)

    map <- setNames(taxon, anno_votu_in$votu_id)
    mat <- aggregate_sum(mat_votu_in, map)

    tax_sum <- anno_votu_in |>
      mutate(taxon = map[votu_id]) |>
      group_by(taxon) |>
      summarise(
        n_member_votu = n(),
        member_votu_ids = collapse_ids(votu_id, MAX_MEMBER_IDS_IN_CELL),
        realm   = safe_mode(realm),
        kingdom = safe_mode(kingdom),
        phylum  = safe_mode(phylum),
        class   = safe_mode(class),
        order   = safe_mode(order),
        family  = safe_mode(family),
        genus   = safe_mode(genus),
        species = safe_mode(species),
        top_host_genus = safe_mode(host_genus),
        top_host_species = safe_mode(host_species),
        mean_cherry_score = if (all(is.na(cherry_score))) NA_real_ else mean(cherry_score, na.rm=TRUE),
        max_cherry_score  = if (all(is.na(cherry_score))) NA_real_ else max(cherry_score, na.rm=TRUE),
        top_method = safe_mode(method),
        .groups="drop"
      ) |>
      mutate(
        feature = taxon,
        level = paste0("ViralTax_", rank),
        display_label = taxon
      ) |>
      select(feature, display_label, level,
             n_member_votu, member_votu_ids,
             realm, kingdom, phylum, class, order, family, genus, species,
             top_host_genus, top_host_species,
             mean_cherry_score, max_cherry_score, top_method)

    map_tbl <- tibble(votu_id = names(map), tax_rank = rank, taxon = unname(map))
    return(list(mat=mat, feat_anno=tax_sum, map_tbl=map_tbl))
  }

  stop("Unknown level: ", level)
}

## ---------- core runner ----------
run_level <- function(level_dirname, level_obj,
                      do_heatmap = TRUE,
                      export_per_sample_all = TRUE,
                      heatmap_label_col = "display_label") {

  level_root <- file.path(OUT_ROOT, level_dirname)
  dir_init(level_root)

  # filter (again at level) with same independent criteria
  f <- filter_matrix(level_obj$mat, meta,
                     min_prev_all = MIN_PREVALENCE_ALL,
                     min_total = MIN_TOTAL_TPM,
                     min_prev_by_group = MIN_PREV_BY_GROUP,
                     min_mean_ra = MIN_MEAN_RA)
  mat_f <- f$mat[, sample_order, drop=FALSE]
  write_tsv0(f$summary, file.path(level_root, "00_manifest", "filter_summary.tsv"))
  write_tsv0(f$diagnostics, file.path(level_root, "00_manifest", "filter_diagnostics.tsv"))
  saveRDS(mat_f, file.path(level_root, "00_manifest", "TPM_matrix_filtered.rds"), compress="xz")

  # align annotations + membership
  feat_anno <- level_obj$feat_anno |> filter(feature %in% rownames(mat_f))
  write_tsv0(feat_anno, file.path(level_root, "00_manifest", "feature_annotation.tsv"))
  write_tsv0(level_obj$map_tbl, file.path(level_root, "00_manifest", "map_membership.tsv"))

  # CLR on TPM
  clr_feat_by_samp <- clr_transform(mat_f, pseudocount = PSEUDOCOUNT)  # features x samples
  saveRDS(clr_feat_by_samp, file.path(level_root, "00_manifest", "CLR_features_by_sample.rds"), compress="xz")

  # RA for per-sample outputs
  ra_feat_by_samp <- sweep(mat_f, 2, colSums(mat_f), "/")
  saveRDS(ra_feat_by_samp, file.path(level_root, "00_manifest", "RA_features_by_sample.rds"), compress="xz")

  for (mname in names(MODELS)) {
    out_dir <- file.path(level_root, mname)
    dir_init(out_dir)

    write_tsv0(tibble(level=level_dirname, model=mname,
                      p_cutoff=P_CUTOFF, fallback_topN=FALLBACK_TOP_N,
                      cluster_rows=CLUSTER_ROWS_DEFAULT, cluster_cols=CLUSTER_COLS_DEFAULT),
               file.path(out_dir, "00_manifest", "params.tsv"))

    # test matrix: samples x features
    Y <- t(clr_feat_by_samp)
    rownames(Y) <- colnames(clr_feat_by_samp)
    colnames(Y) <- rownames(clr_feat_by_samp)

    if (MODELS[[mname]]$type == "resid") {
      Y <- residualize_matrix(Y, MODELS[[mname]]$X)
    }

    # per-feature Wilcoxon
    idx_gr <- which(meta$group=="PD_GR")
    idx_mr <- which(meta$group=="PD_MR")
    feats <- colnames(Y)

    st <- t(vapply(feats, function(f) wilcox_one(Y[idx_mr, f], Y[idx_gr, f]), numeric(5)))

    res <- tibble(
      feature = feats,
      effect_median_diff = st[, "median_diff"],
      W = st[, "W"],
      U = st[, "U"],
      rank_biserial = st[, "rank_biserial"],
      p_value = st[, "p_value"]
    ) |>
      mutate(q_value = p.adjust(p_value, method="BH"),
             direction = case_when(effect_median_diff>0 ~ "PD_MR>PD_GR",
                                   effect_median_diff<0 ~ "PD_MR<PD_GR",
                                   TRUE ~ "NA")) |>
      left_join(feat_anno, by="feature")

    write_tsv0(res, file.path(out_dir, "tables", "all_results.tsv"))
    write_tsv0(res |> filter(q_value < FDR_CUTOFF),
               file.path(out_dir, "tables", paste0("significant_FDR_",FDR_CUTOFF,".tsv")))

    # selection for heatmap (p<0.05 or fallback top N)
    sel <- res |> filter(p_value < P_CUTOFF) |> arrange(p_value)
    note <- ""
    if (nrow(sel) == 0) {
      sel <- res |> arrange(p_value) |> slice_head(n = min(FALLBACK_TOP_N, nrow(res)))
      note <- paste0("No features with p<", P_CUTOFF, "; fallback top ", nrow(sel), " by p.")
    } else {
      note <- paste0("Selected features: p<", P_CUTOFF, " (n=", nrow(sel), ").")
    }
    writeLines(note, con = file.path(out_dir, "logs", "selection_note.txt"))
    write_tsv0(sel, file.path(out_dir, "tables", "selected_by_p.tsv"))

    # test_value matrix (features x samples)
    test_feat_by_samp <- t(Y)
    test_feat_by_samp <- test_feat_by_samp[, sample_order, drop=FALSE]
    # ---- Host wide tables requested (feature rows, sample columns + stats) ----
    if (level_dirname %in% c("Host", "Host_species")) {

      prefix <- if (level_dirname == "Host") "HostGenus" else "HostSpecies"

      out_all <- file.path(out_dir, "tables", paste0(prefix, "_testvalue_wide_with_stats.tsv"))
      out_p05 <- file.path(out_dir, "tables", paste0(prefix, "_testvalue_wide_with_stats_p_lt_0.05.tsv"))

      write_feature_sample_wide_with_stats(
        test_feat_by_samp = test_feat_by_samp,
        res_df = res,
        sample_order = sample_order,
        out_all_path = out_all,
        out_p05_path = out_p05,
        p_cutoff = 0.05
      )
    }

    # ---- per-sample ALL (WIDE as requested) ----
    # The file name per_sample_all_long.tsv is kept, but content is WIDE (sample rows).
    if (export_per_sample_all) {
      # wide matrices
      write_matrix_wide_sample_rows(mat_f, file.path(out_dir, "tables", "per_sample_all_TPM_wide.tsv"), sample_order)
      write_matrix_wide_sample_rows(ra_feat_by_samp, file.path(out_dir, "tables", "per_sample_all_RA_wide.tsv"), sample_order)
      write_matrix_wide_sample_rows(clr_feat_by_samp, file.path(out_dir, "tables", "per_sample_all_CLR_wide.tsv"), sample_order)
      write_matrix_wide_sample_rows(test_feat_by_samp, file.path(out_dir, "tables", "per_sample_all_testvalue_wide.tsv"), sample_order)

      # metadata + testvalue combined wide (main “single-sample table”)
      # --- per_sample_all_long.tsv as WIDE (feature rows, sample columns) ---
      tv_feat <- as.data.frame(test_feat_by_samp) |> tibble::rownames_to_column("feature")
      # join annotations (keep as many annotation columns as possible)
      tv_feat_anno <- feat_anno |> distinct(feature, .keep_all = TRUE) |> left_join(tv_feat, by = "feature")

      # reorder: annotation columns first, then sample columns
      ann_cols <- setdiff(names(tv_feat_anno), sample_order)
      tv_feat_anno <- tv_feat_anno[, c(ann_cols, sample_order), drop = FALSE]

      write_tsv0(tv_feat_anno, file.path(out_dir, "tables", "per_sample_all_long.tsv"))

      # also export tidy long for completeness (optional)
      tv_long <- as.data.frame(test_feat_by_samp) |> rownames_to_column("feature") |>
        pivot_longer(-feature, names_to="sample_id", values_to="test_value") |>
        left_join(meta, by="sample_id")
      write_tsv0(tv_long, file.path(out_dir, "tables", "per_sample_all_long_tidy.tsv"))
    }

    # ---- per-sample SELECTED (both long + heatmap matrices) ----
    sel_feats <- sel$feature
    tv_sel  <- test_feat_by_samp[sel_feats, , drop=FALSE]
    z_sel   <- row_zscore(tv_sel)

    write_tsv0(as.data.frame(tv_sel) |> rownames_to_column("feature"),
               file.path(out_dir, "tables", "heatmap_matrix_testvalue_wide.tsv"))
    write_tsv0(as.data.frame(z_sel) |> rownames_to_column("feature"),
               file.path(out_dir, "tables", "heatmap_matrix_zscore_wide.tsv"))
    write_tsv0(meta, file.path(out_dir, "tables", "sample_annotation.tsv"))

    # selected long (for plotting in any software)
    sel_long <- as.data.frame(tv_sel) |> rownames_to_column("feature") |>
      pivot_longer(-feature, names_to="sample_id", values_to="test_value") |>
      left_join(as.data.frame(z_sel) |> rownames_to_column("feature") |>
                  pivot_longer(-feature, names_to="sample_id", values_to="z_value"),
                by=c("feature","sample_id")) |>
      left_join(meta, by="sample_id") |>
      left_join(sel |> select(feature, p_value, q_value, effect_median_diff, rank_biserial, direction,
                              all_of(intersect(names(sel), c("display_label","level","top_host_genus","top_host_species",
                                                             "realm","kingdom","phylum","class","order","family","genus","species",
                                                             "top_phylum","top_family","top_host_species","top_method","n_member_votu")))),
                by="feature")
    write_tsv0(sel_long, file.path(out_dir, "tables", "per_sample_selected_long.tsv"))

    # ---- heatmap plot (PDF required) ----
    if (do_heatmap) {
      # display labels for rows if available
      rowlabs <- sel_feats
      if (!is.null(heatmap_label_col) && heatmap_label_col %in% names(sel)) {
        rowlabs <- sel[[heatmap_label_col]]
        if (any(duplicated(rowlabs))) {
          rowlabs <- paste0(rowlabs, " [", sel_feats, "]")
        }
      }
      mat_plot <- z_sel
      rownames(mat_plot) <- rowlabs

      anno_col <- data.frame(group = meta$group); rownames(anno_col) <- meta$sample_id
      anno_col <- anno_col[sample_order, , drop=FALSE]

      out_pdf <- file.path(out_dir, "figures", paste0("heatmap_", level_dirname, "_", mname, ".pdf"))
      out_png <- file.path(out_dir, "figures", paste0("heatmap_", level_dirname, "_", mname, ".png"))
      plot_heatmap(mat_plot, anno_col, out_pdf, out_png,
                   cluster_rows = CLUSTER_ROWS_DEFAULT,
                   cluster_cols = CLUSTER_COLS_DEFAULT)
    }

    capture.output(sessionInfo(), file = file.path(out_dir, "logs", "sessionInfo.txt"))
    message("Done: ", level_dirname, " | ", mname, " | n=", nrow(res), " | heatmap=", do_heatmap,
            " | per_sample_all_wide=", export_per_sample_all)
  }
}

## ---------- RUN ----------
# vOTU stats only (no heatmap)
if (DO_vOTU_LEVEL) {
  lvl_votu <- build_level("vOTU", mat_votu_in = mat_votu0, anno_votu_in = anno_votu0)
  run_level("vOTU", lvl_votu,
            do_heatmap = FALSE,
            export_per_sample_all = vOTU_EXPORT_PER_SAMPLE_ALL,
            heatmap_label_col = "display_label")
}

# Host genus + Host species
lvl_host_genus <- build_level("Host", mat_votu_in = mat_votu0, anno_votu_in = anno_votu0)
run_level("Host", lvl_host_genus,
          do_heatmap = MAKE_HEATMAP_HOST,
          export_per_sample_all = TRUE,
          heatmap_label_col = "display_label")

lvl_host_species <- build_level("Host_species", mat_votu_in = mat_votu0, anno_votu_in = anno_votu0)
run_level("Host_species", lvl_host_species,
          do_heatmap = MAKE_HEATMAP_HOST,
          export_per_sample_all = TRUE,
          heatmap_label_col = "display_label")

# Viral taxonomy ranks
for (rk in VIRAL_RANKS) {
  lvl_tax <- build_level("ViralTax", rank = rk, mat_votu_in = mat_votu0, anno_votu_in = anno_votu0)
  run_level(paste0("ViralTax_", rk), lvl_tax,
            do_heatmap = MAKE_HEATMAP_TAX,
            export_per_sample_all = TRUE,
            heatmap_label_col = "display_label")
}

## ---------- vOTU presence supplement ----------
if (RUN_vOTU_PRESENCE_SUPPLEMENT) {
  out_dir <- file.path(OUT_ROOT, "vOTU", "presence_supplement")
  dir_init(out_dir)

  prev_all <- rowMeans(mat_votu > 0)
  keep <- prev_all >= PRES_MIN_PREVALENCE_ALL
  mat_p <- mat_votu[keep, sample_order, drop=FALSE]

  write_tsv0(tibble(n_total=nrow(mat_votu), n_keep=nrow(mat_p),
                    pres_min_prev_all=PRES_MIN_PREVALENCE_ALL,
                    pres_test_method=PRES_TEST_METHOD),
             file.path(out_dir, "00_manifest", "presence_params.tsv"))

  if (nrow(mat_p) > 0) {
    idx_gr <- which(meta$group=="PD_GR")
    idx_mr <- which(meta$group=="PD_MR")
    n_gr <- length(idx_gr); n_mr <- length(idx_mr)

    pres <- (mat_p > 0) * 1
    a <- rowSums(pres[, idx_mr, drop=FALSE])
    c <- rowSums(pres[, idx_gr, drop=FALSE])
    b <- n_mr - a
    d <- n_gr - c

    OR <- ((a+0.5)*(d+0.5))/((b+0.5)*(c+0.5))
    logOR <- log(OR)

    method_use <- PRES_TEST_METHOD
    if (method_use=="auto") method_use <- if (nrow(mat_p) <= 5000) "fisher" else "chisq"

    if (method_use=="chisq") {
      n <- n_mr + n_gr
      denom <- (a+b)*(c+d)*(a+c)*(b+d)
      chi2 <- (n*(a*d - b*c)^2)/pmax(denom, 1e-12)
      p <- pchisq(chi2, df=1, lower.tail=FALSE)
    } else {
      p <- vapply(seq_len(nrow(mat_p)), function(i) {
        m <- matrix(c(a[i], b[i], c[i], d[i]), nrow=2, byrow=TRUE)
        suppressWarnings(fisher.test(m)$p.value)
      }, numeric(1))
    }

    q <- p.adjust(p, method="BH")

    pres_res <- tibble(
      feature = rownames(mat_p),
      MR_present=a, MR_absent=b, GR_present=c, GR_absent=d,
      prev_MR=a/n_mr, prev_GR=c/n_gr,
      logOR=logOR, OR=OR,
      p_value=p, q_value=q,
      direction = case_when(logOR>0 ~ "PD_MR>PD_GR",
                            logOR<0 ~ "PD_MR<PD_GR",
                            TRUE ~ "NA")
    ) |>
      left_join(anno_votu |> rename(feature=votu_id), by="feature")

    write_tsv0(pres_res, file.path(out_dir, "tables", "vOTU_presence_all_results.tsv"))
    write_tsv0(pres_res |> filter(q_value < FDR_CUTOFF),
               file.path(out_dir, "tables", paste0("vOTU_presence_significant_FDR_",FDR_CUTOFF,".tsv")))
  }
}

message("All virome DA + heatmaps done -> ", OUT_ROOT)