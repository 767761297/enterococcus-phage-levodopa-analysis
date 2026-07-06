## STEP2_heatmap_signature_q0.1_6levels.R
## 6 taxonomic levels fingerprint heatmaps (4 groups in one figure)
## Input:
##   1) Signature tables from Step1: 16S_Signature_q0.1_6Levels/SIGNATURE_q0.1_<level>_GR_vs_MR.csv
##   2) Abundance tables: allsample_bundle_from_featuretable/allsample.<level>.abundance.xls
## Output:
##   - 6 PNG + 6 PDF (one per level)
##   - one merged multipage PDF

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(tibble)
  library(pheatmap)
  library(grid)
})

# -------------------- SETTINGS --------------------
SIG_DIR <- "4. 16S_Signature_q0.1_6Levels"
IN_DIR  <- "allsample_bundle_from_featuretable"
OUTDIR  <- "16S_FingerprintHeatmaps_q0.1_6Levels"
dir.create(OUTDIR, showWarnings = FALSE)

LEVELS <- c("phylum","class","order","family","genus","species")

# Group mapping
map_group <- function(s) {
  if (str_detect(s, "^A")) return("GR")
  if (str_detect(s, "^B")) return("MR")
  if (str_detect(s, "^C")) return("MR+Phage cocktail")
  if (str_detect(s, "^D")) return("MR+AFMT")
  return(NA_character_)
}
GROUP_LEVELS <- c("GR","MR","MR+Phage cocktail","MR+AFMT")

# Column order in heatmap (recommended for interpretation)
SAMPLE_ORDER <- c(paste0("A",1:5), paste0("B",1:5), paste0("C",1:5), paste0("D",1:5))

# Use CLR on counts (counts + pseudo), then heatmap shows row-zscore
PSEUDO <- 1

# Heatmap options
CLUSTER_COLS <- FALSE      # keep group order
SHOW_COLNAMES <- FALSE     # usually cleaner
FONTSIZE_ROW <- 7
# --------------------------------------------------

read_abundance <- function(level) {
  f <- file.path(IN_DIR, paste0("allsample.", level, ".abundance.xls"))
  if (!file.exists(f)) stop("Missing abundance file: ", f)

  df <- read.delim(f, sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
  tax_col <- colnames(df)[1]
  df <- df %>% rename(taxon = all_of(tax_col))

  samp <- intersect(colnames(df), SAMPLE_ORDER)
  if (length(samp) == 0) {
    # fallback: auto-detect A1/B2/C3/D4 columns
    samp <- setdiff(colnames(df), "taxon")
    samp <- samp[grepl("^[A-D][0-9]+$", samp)]
  }
  if (length(samp) == 0) stop("No sample columns detected in: ", f)

  mat <- as.matrix(df[, samp, drop=FALSE])
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  rownames(mat) <- make.unique(as.character(df$taxon))

  # reorder columns to SAMPLE_ORDER (keep those that exist)
  samp2 <- SAMPLE_ORDER[SAMPLE_ORDER %in% colnames(mat)]
  mat <- mat[, samp2, drop=FALSE]

  mat
}

clr_transform <- function(counts, pseudo=1) {
  x <- counts + pseudo
  logx <- log(x)
  sweep(logx, 2, colMeans(logx), "-")
}

read_signature_taxa <- function(level) {
  f <- file.path(SIG_DIR, paste0("SIGNATURE_q0.1_", level, "_GR_vs_MR.csv"))
  if (!file.exists(f)) stop("Missing signature table: ", f)
  sig <- readr::read_csv(f, show_col_types = FALSE)
  if (!("taxon" %in% names(sig))) stop("Signature file lacks 'taxon' column: ", f)
  unique(as.character(sig$taxon))
}

make_heatmap_one_level <- function(level) {
  sig_taxa <- read_signature_taxa(level)
  if (length(sig_taxa) == 0) {
    message("[SKIP] ", level, ": signature set is empty (q<0.1 produced 0 taxa).")
    return(invisible(NULL))
  }

  counts <- read_abundance(level)

  keep <- intersect(sig_taxa, rownames(counts))
  if (length(keep) == 0) {
    message("[SKIP] ", level, ": none of signature taxa found in abundance table.")
    return(invisible(NULL))
  }

  counts_sig <- counts[keep, , drop=FALSE]
  clr_sig <- clr_transform(counts_sig, PSEUDO)

  # order columns + meta
  samples <- colnames(clr_sig)
  groups  <- vapply(samples, map_group, character(1))
  groups  <- factor(groups, levels = GROUP_LEVELS)

  # annotation + gaps (only meaningful if not clustering cols)
  anno_col <- data.frame(group = groups)
  rownames(anno_col) <- samples
  gaps_col <- cumsum(as.integer(table(groups)))

  # remove rows with zero variance (avoid NA z-score)
  row_sd <- apply(clr_sig, 1, sd)
  clr_sig <- clr_sig[row_sd > 0, , drop=FALSE]
  if (nrow(clr_sig) == 0) {
    message("[SKIP] ", level, ": all signature taxa have zero variance after filtering.")
    return(invisible(NULL))
  }

  # dynamic figure height
  h_in <- max(4, min(14, 0.22 * nrow(clr_sig) + 2.5))
  w_in <- 9.5

  out_png <- file.path(OUTDIR, paste0("Heatmap_signature_q0.1_", level, "_4groups.png"))
  out_pdf <- file.path(OUTDIR, paste0("Heatmap_signature_q0.1_", level, "_4groups.pdf"))

  # IMPORTANT: main must NOT be NULL (pheatmap 会报你之前遇到的 length=0 错误)
  # Use main=NA to remove title cleanly
  pheatmap(
    clr_sig,
    scale = "row",
    annotation_col = anno_col,
    cluster_rows = TRUE,
    cluster_cols = CLUSTER_COLS,
    gaps_col = if (!CLUSTER_COLS) gaps_col else NULL,
    show_colnames = SHOW_COLNAMES,
    fontsize_row = FONTSIZE_ROW,
    fontsize_col = 8,
    main = NA,
    filename = out_png,
    width = w_in, height = h_in
  )

  pheatmap(
    clr_sig,
    scale = "row",
    annotation_col = anno_col,
    cluster_rows = TRUE,
    cluster_cols = CLUSTER_COLS,
    gaps_col = if (!CLUSTER_COLS) gaps_col else NULL,
    show_colnames = SHOW_COLNAMES,
    fontsize_row = FONTSIZE_ROW,
    fontsize_col = 8,
    main = NA,
    filename = out_pdf,
    width = w_in, height = h_in
  )

  message(sprintf("[OK] %s: taxa=%d  -> %s", level, nrow(clr_sig), out_png))
  invisible(list(level=level, n_taxa=nrow(clr_sig), png=out_png, pdf=out_pdf))
}

# ---- run 6 levels ----
res <- lapply(LEVELS, make_heatmap_one_level)

# ---- merged multipage PDF (one page per level) ----
merged_pdf <- file.path(OUTDIR, "Heatmap_signature_q0.1_6levels_merged.pdf")
pdf(merged_pdf, width = 10, height = 8)
for (lv in LEVELS) {
  sig_taxa <- tryCatch(read_signature_taxa(lv), error=function(e) character(0))
  if (length(sig_taxa) == 0) next

  counts <- tryCatch(read_abundance(lv), error=function(e) NULL)
  if (is.null(counts)) next
  keep <- intersect(sig_taxa, rownames(counts))
  if (length(keep) == 0) next

  clr_sig <- clr_transform(counts[keep, , drop=FALSE], PSEUDO)
  row_sd <- apply(clr_sig, 1, sd)
  clr_sig <- clr_sig[row_sd > 0, , drop=FALSE]
  if (nrow(clr_sig) == 0) next

  samples <- colnames(clr_sig)
  groups  <- factor(vapply(samples, map_group, character(1)), levels = GROUP_LEVELS)
  anno_col <- data.frame(group = groups)
  rownames(anno_col) <- samples
  gaps_col <- cumsum(as.integer(table(groups)))

  ph <- pheatmap(
    clr_sig,
    scale = "row",
    annotation_col = anno_col,
    cluster_rows = TRUE,
    cluster_cols = CLUSTER_COLS,
    gaps_col = if (!CLUSTER_COLS) gaps_col else NULL,
    show_colnames = SHOW_COLNAMES,
    fontsize_row = FONTSIZE_ROW,
    fontsize_col = 8,
    main = NA
  )
  grid.newpage()
  grid.draw(ph$gtable)
}
dev.off()

message("\n[DONE] Heatmaps saved in: ", OUTDIR)
message("[DONE] Merged PDF: ", merged_pdf)
