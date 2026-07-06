## ============================================================
## STEP_16S_stats_from_allsample_bundle_6levels_vsMR_q0.1.R
## Input: allsample_bundle_from_featuretable/allsample.<rank>.reabundance.xls
## Output: one Excel per rank:
##   - all_results (full table)
##   - signature_q0.1 (q<0.1 in any vsMR comparison)
##   - sig_GR_vs_MR / sig_MR_Phage_vs_MR / sig_MR_AFMT_vs_MR
##   - meta
##
## Stats: CLR + Wilcoxon (vs MR) + BH(FDR), plus Kruskal-Wallis overall p/fdr
## log2FC: log2(geom-mean(rel_abund + pseudo) ratio) vs MR (platform-like)
## ============================================================

options(stringsAsFactors = FALSE)

## -----------------------------
## 0) USER PARAMETERS (edit here)
## -----------------------------
BUNDLE_DIR <- "allsample_bundle_from_featuretable"   # folder in your working directory
OUTDIR     <- "16S_STATS_vsMR_q0.1_6levels_fromBundle"

Q_CUTOFF   <- 0.10
PSEUDO_REL <- 1e-6  # pseudo-count for rel abundance (CLR + geom mean)

## Sample group mapping by sample ID prefix
## Default: A*=GR, B*=MR, C*=MR+Phage, D*=MR+AFMT
GROUP_PATTERNS <- list(
  GR       = "^A",
  MR       = "^B",
  MR_Phage = "^C",
  MR_AFMT  = "^D"
)

GROUP_DISPLAY <- c(
  GR       = "GR",
  MR       = "MR",
  MR_Phage = "MR+Phage cocktail",
  MR_AFMT  = "MR+AFMT"
)

## -----------------------------
## 1) Dependencies
## -----------------------------
if (!requireNamespace("writexl", quietly = TRUE)) {
  install.packages("writexl", repos = "https://cloud.r-project.org")
}

## -----------------------------
## 2) Helper functions
## -----------------------------
trim <- function(x) gsub("^\\s+|\\s+$", "", x)

safe_make_unique_names <- function(x) {
  x <- trim(x)
  x[x == ""] <- NA_character_
  x[is.na(x)] <- paste0("X", seq_len(sum(is.na(x))))
  make.unique(x)
}

make_group_vector <- function(sample_names, patterns_named_list) {
  grp <- rep(NA_character_, length(sample_names))
  for (g in names(patterns_named_list)) {
    pat <- patterns_named_list[[g]]
    hit <- grepl(pat, sample_names, perl = TRUE)
    grp[hit] <- g
  }
  grp
}

clr_transform <- function(rel_mat, pseudo = 1e-6) {
  x  <- rel_mat + pseudo
  lx <- log(x)
  gm <- colMeans(lx, na.rm = TRUE)
  sweep(lx, 2, gm, "-")
}

geom_mean <- function(x, pseudo = 1e-6) {
  x <- x + pseudo
  exp(mean(log(x), na.rm = TRUE))
}

safe_wilcox_p <- function(x, y) {
  tryCatch(wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA_real_)
}

safe_kw_p <- function(x, g) {
  tryCatch(kruskal.test(x ~ g)$p.value, error = function(e) NA_real_)
}

read_allsample_reabundance <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  df <- read.delim(path, header = TRUE, sep = "\t", quote = "", comment.char = "",
                   check.names = FALSE, stringsAsFactors = FALSE)
  names(df) <- safe_make_unique_names(names(df))
  if (ncol(df) < 3) stop("Unexpected table format (need taxon + >=2 samples): ", path)

  tax_col <- names(df)[1]
  df[[tax_col]] <- as.character(df[[tax_col]])
  names(df)[1] <- "taxon"

  ## Coerce sample columns to numeric
  for (j in 2:ncol(df)) df[[j]] <- suppressWarnings(as.numeric(df[[j]]))
  df
}

## -----------------------------
## 3) Main per-rank function
## -----------------------------
run_rank <- function(rank) {

  in_file <- file.path(BUNDLE_DIR, paste0("allsample.", rank, ".reabundance.xls"))
  df <- read_allsample_reabundance(in_file)

  taxa <- df$taxon
  sample_names <- names(df)[-1]

  mat_rel <- as.matrix(df[, -1, drop = FALSE])
  rownames(mat_rel) <- taxa
  colnames(mat_rel) <- sample_names

  ## group vector
  group_vec <- make_group_vector(sample_names, GROUP_PATTERNS)
  if (any(is.na(group_vec))) {
    bad <- sample_names[is.na(group_vec)]
    stop("Cannot assign group for these sample columns: ", paste(bad, collapse = ", "),
         "\nPlease edit GROUP_PATTERNS.")
  }

  ## Ensure MR as reference exists
  if (sum(group_vec == "MR") == 0) stop("No MR samples detected. Check GROUP_PATTERNS.")

  group_vec <- factor(group_vec, levels = c("MR","GR","MR_Phage","MR_AFMT"))

  idx_MR    <- which(group_vec == "MR")
  idx_GR    <- which(group_vec == "GR")
  idx_Phage <- which(group_vec == "MR_Phage")
  idx_AFMT  <- which(group_vec == "MR_AFMT")

  ## CLR
  mat_clr <- clr_transform(mat_rel, pseudo = PSEUDO_REL)

  ## group means/sd (relative abundance)
  mean_MR    <- rowMeans(mat_rel[, idx_MR, drop=FALSE], na.rm = TRUE)
  sd_MR      <- apply(mat_rel[, idx_MR, drop=FALSE], 1, sd, na.rm = TRUE)

  mean_GR    <- rowMeans(mat_rel[, idx_GR, drop=FALSE], na.rm = TRUE)
  sd_GR      <- apply(mat_rel[, idx_GR, drop=FALSE], 1, sd, na.rm = TRUE)

  mean_Phage <- rowMeans(mat_rel[, idx_Phage, drop=FALSE], na.rm = TRUE)
  sd_Phage   <- apply(mat_rel[, idx_Phage, drop=FALSE], 1, sd, na.rm = TRUE)

  mean_AFMT  <- rowMeans(mat_rel[, idx_AFMT, drop=FALSE], na.rm = TRUE)
  sd_AFMT    <- apply(mat_rel[, idx_AFMT, drop=FALSE], 1, sd, na.rm = TRUE)

  ## geometric means (relative abundance)
  gm_MR    <- apply(mat_rel[, idx_MR, drop=FALSE], 1, geom_mean, pseudo = PSEUDO_REL)
  gm_GR    <- apply(mat_rel[, idx_GR, drop=FALSE], 1, geom_mean, pseudo = PSEUDO_REL)
  gm_Phage <- apply(mat_rel[, idx_Phage, drop=FALSE], 1, geom_mean, pseudo = PSEUDO_REL)
  gm_AFMT  <- apply(mat_rel[, idx_AFMT, drop=FALSE], 1, geom_mean, pseudo = PSEUDO_REL)

  ## log2FC vs MR (gm ratio)
  log2FC_GR_vs_MR    <- log2(gm_GR / gm_MR)
  log2FC_Phage_vs_MR <- log2(gm_Phage / gm_MR)
  log2FC_AFMT_vs_MR  <- log2(gm_AFMT / gm_MR)

  ## p-values: Wilcoxon on CLR
  p_GR_vs_MR <- vapply(seq_len(nrow(mat_clr)), function(i)
    safe_wilcox_p(mat_clr[i, idx_GR], mat_clr[i, idx_MR]), numeric(1))

  p_Phage_vs_MR <- vapply(seq_len(nrow(mat_clr)), function(i)
    safe_wilcox_p(mat_clr[i, idx_Phage], mat_clr[i, idx_MR]), numeric(1))

  p_AFMT_vs_MR <- vapply(seq_len(nrow(mat_clr)), function(i)
    safe_wilcox_p(mat_clr[i, idx_AFMT], mat_clr[i, idx_MR]), numeric(1))

  ## overall p: Kruskal-Wallis on CLR
  p_overall <- vapply(seq_len(nrow(mat_clr)), function(i)
    safe_kw_p(as.numeric(mat_clr[i, ]), group_vec), numeric(1))

  ## BH FDR within each contrast
  fdr_GR_vs_MR    <- p.adjust(p_GR_vs_MR,    method = "BH")
  fdr_Phage_vs_MR <- p.adjust(p_Phage_vs_MR, method = "BH")
  fdr_AFMT_vs_MR  <- p.adjust(p_AFMT_vs_MR,  method = "BH")
  fdr_overall     <- p.adjust(p_overall,     method = "BH")

  ## full table (metabolomics-like)
  all_results <- data.frame(
    taxon = taxa,
    p_overall = p_overall,

    mean_GR = mean_GR, sd_GR = sd_GR,
    mean_MR = mean_MR, sd_MR = sd_MR,
    mean_MR_Phage = mean_Phage, sd_MR_Phage = sd_Phage,
    mean_MR_AFMT  = mean_AFMT,  sd_MR_AFMT  = sd_AFMT,

    log2FC_GR_vs_MR = log2FC_GR_vs_MR,
    p_GR_vs_MR      = p_GR_vs_MR,

    log2FC_MR_Phage_vs_MR = log2FC_Phage_vs_MR,
    p_MR_Phage_vs_MR      = p_Phage_vs_MR,

    log2FC_MR_AFMT_vs_MR = log2FC_AFMT_vs_MR,
    p_MR_AFMT_vs_MR      = p_AFMT_vs_MR,

    fdr_GR_vs_MR        = fdr_GR_vs_MR,
    fdr_MR_Phage_vs_MR  = fdr_Phage_vs_MR,
    fdr_MR_AFMT_vs_MR   = fdr_AFMT_vs_MR,
    fdr_overall         = fdr_overall,

    gm_GR = gm_GR,
    gm_MR = gm_MR,
    gm_MR_Phage = gm_Phage,
    gm_MR_AFMT  = gm_AFMT,
    stringsAsFactors = FALSE
  )

  ## signature: q<0.1 in any of the 3 vsMR contrasts
  sig_GR    <- all_results$fdr_GR_vs_MR       < Q_CUTOFF
  sig_Phage <- all_results$fdr_MR_Phage_vs_MR < Q_CUTOFF
  sig_AFMT  <- all_results$fdr_MR_AFMT_vs_MR  < Q_CUTOFF

  signature_q01 <- within(all_results, {
    sig_GR_vs_MR       <- sig_GR
    sig_MR_Phage_vs_MR <- sig_Phage
    sig_MR_AFMT_vs_MR  <- sig_AFMT
    sig_any            <- sig_GR_vs_MR | sig_MR_Phage_vs_MR | sig_MR_AFMT_vs_MR

    dir_GR_vs_MR       <- ifelse(log2FC_GR_vs_MR > 0, "UP", "DOWN")
    dir_MR_Phage_vs_MR <- ifelse(log2FC_MR_Phage_vs_MR > 0, "UP", "DOWN")
    dir_MR_AFMT_vs_MR  <- ifelse(log2FC_MR_AFMT_vs_MR > 0, "UP", "DOWN")
  })
  signature_q01 <- signature_q01[signature_q01$sig_any, , drop = FALSE]

  ## sorted for readability
  signature_q01 <- signature_q01[order(signature_q01$fdr_GR_vs_MR,
                                       signature_q01$fdr_MR_Phage_vs_MR,
                                       signature_q01$fdr_MR_AFMT_vs_MR,
                                       signature_q01$taxon), , drop = FALSE]

  sig_GR_tbl    <- all_results[sig_GR, , drop = FALSE]
  sig_Phage_tbl <- all_results[sig_Phage, , drop = FALSE]
  sig_AFMT_tbl  <- all_results[sig_AFMT, , drop = FALSE]

  meta <- data.frame(
    key = c("rank","input_file","pseudo_rel","q_cutoff",
            "pairwise_test","overall_test","p_adjust",
            "groups",
            "n_GR","n_MR","n_MR_Phage","n_MR_AFMT"),
    value = c(rank, in_file, as.character(PSEUDO_REL), as.character(Q_CUTOFF),
              "Wilcoxon on CLR", "Kruskal-Wallis on CLR", "BH (FDR)",
              paste(names(GROUP_DISPLAY), "=", GROUP_DISPLAY[names(GROUP_DISPLAY)], collapse = "; "),
              sum(group_vec=="GR"), sum(group_vec=="MR"),
              sum(group_vec=="MR_Phage"), sum(group_vec=="MR_AFMT")),
    stringsAsFactors = FALSE
  )

  list(
    meta = meta,
    all_results = all_results,
    `signature_q0.1` = signature_q01,
    sig_GR_vs_MR = sig_GR_tbl,
    sig_MR_Phage_vs_MR = sig_Phage_tbl,
    sig_MR_AFMT_vs_MR = sig_AFMT_tbl
  )
}

## -----------------------------
## 4) Run 6 ranks & write Excels
## -----------------------------
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

ranks <- c("phylum","class","order","family","genus","species")

for (rk in ranks) {
  message("[RUN] ", rk)
  wb <- run_rank(rk)
  out_xlsx <- file.path(OUTDIR, paste0("16S_", rk, "_stats_vsMR_q0.1.xlsx"))
  writexl::write_xlsx(wb, out_xlsx)
  message("  [OK] wrote: ", out_xlsx,
          " | taxa=", nrow(wb$all_results),
          " | signature_any(q<0.1)=", nrow(wb$`signature_q0.1`))
}

message("[DONE] All 6 workbooks saved in: ", OUTDIR)
