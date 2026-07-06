## STEP1_build_16S_signature_q0.1_6levels.R
## Step 1 only: build fingerprint/signature sets (q < 0.1) at 6 taxonomic levels
## Method: CLR (counts + 1) + Wilcoxon (GR vs MR) + BH(FDR)
## Input: allsample.<level>.abundance.xls under IN_DIR (TSV .xls is OK)

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(tibble)
})

# -------------------- USER SETTINGS --------------------
IN_DIR <- "allsample_bundle_from_featuretable"   # where allsample.*.abundance.xls are
OUTDIR <- "16S_Signature_q0.1_6Levels"
dir.create(OUTDIR, showWarnings = FALSE)

LEVELS <- c("phylum","class","order","family","genus","species")

# your sample IDs
SAMPLE_ORDER <- c(paste0("A",1:5), paste0("C",1:5), paste0("D",1:5), paste0("B",1:5))

# group mapping (only GR/MR used here)
map_group <- function(s) {
  if (str_detect(s, "^A")) return("GR")
  if (str_detect(s, "^B")) return("MR")
  if (str_detect(s, "^C")) return("MR+Phage cocktail")
  if (str_detect(s, "^D")) return("MR+AFMT")
  return(NA_character_)
}

# Stats settings
PSEUDO <- 1          # for CLR on counts
PREV_CUT <- 0.10     # prevalence filter across ALL samples (count>0 fraction). set 0 to disable
Q_CUT <- 0.10        # signature threshold
# ------------------------------------------------------


# -------------------- HELPERS --------------------
read_abundance_table <- function(level) {
  f1 <- file.path(IN_DIR, paste0("allsample.", level, ".abundance.xls"))
  f2 <- paste0("allsample.", level, ".abundance.xls")  # fallback: current dir

  f <- if (file.exists(f1)) f1 else f2
  if (!file.exists(f)) stop("Cannot find input file for level '", level, "': ", f1, " or ", f2)

  df <- read.delim(f, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(ncol(df) >= 2)

  tax_col <- colnames(df)[1]
  df <- df %>% rename(taxon = all_of(tax_col))

  # pick sample columns present in table
  sample_cols <- intersect(colnames(df), SAMPLE_ORDER)
  if (length(sample_cols) == 0) {
    # fallback: detect A1/B2/C3/D4-like columns
    sample_cols <- setdiff(colnames(df), "taxon")
    sample_cols <- sample_cols[grepl("^[A-D][0-9]+$", sample_cols)]
  }
  if (length(sample_cols) == 0) stop("No sample columns detected in: ", f)

  # numeric matrix
  mat <- as.matrix(df[, sample_cols, drop = FALSE])
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0

  rownames(mat) <- make.unique(as.character(df$taxon))
  list(mat = mat, file = f)
}

clr_transform <- function(counts, pseudo = 1) {
  x <- counts + pseudo
  logx <- log(x)
  sweep(logx, 2, colMeans(logx), "-")
}

safe_wilcox <- function(x1, x0) {
  if (length(unique(c(x1, x0))) <= 1) return(1)
  out <- tryCatch(
    suppressWarnings(wilcox.test(x1, x0, exact = FALSE)$p.value),
    error = function(e) NA_real_
  )
  out
}

calc_log2fc_rel <- function(counts, s1, s0, eps = 1e-6) {
  rel <- sweep(counts, 2, colSums(counts), "/")
  m1 <- rowMeans(rel[, s1, drop = FALSE]) + eps
  m0 <- rowMeans(rel[, s0, drop = FALSE]) + eps
  log2(m1 / m0)
}

# -------------------- MAIN --------------------
summary_rows <- list()

for (lv in LEVELS) {
  message("\n[RUN] level = ", lv)

  obj <- read_abundance_table(lv)
  counts <- obj$mat

  # meta
  samples <- colnames(counts)
  meta <- tibble(sample = samples, group = vapply(samples, map_group, character(1))) %>%
    filter(group %in% c("GR", "MR", "MR+Phage cocktail", "MR+AFMT"))

  # Keep only samples we can map (should be all)
  keep_samp <- meta$sample
  counts <- counts[, keep_samp, drop = FALSE]

  # prevalence filter across all samples
  if (PREV_CUT > 0) {
    prv <- rowMeans(counts > 0)
    counts <- counts[prv >= PREV_CUT, , drop = FALSE]
  }

  # GR vs MR only
  s_GR <- meta$sample[meta$group == "GR"]
  s_MR <- meta$sample[meta$group == "MR"]
  if (length(s_GR) < 2 || length(s_MR) < 2) stop("Not enough samples in GR or MR for: ", lv)

  clr <- clr_transform(counts, pseudo = PSEUDO)

  delta_clr <- rowMeans(clr[, s_GR, drop = FALSE]) - rowMeans(clr[, s_MR, drop = FALSE])
  pvals <- apply(clr, 1, function(v) safe_wilcox(v[s_GR], v[s_MR]))
  qvals <- p.adjust(pvals, method = "BH")
  log2fc_rel <- calc_log2fc_rel(counts, s_GR, s_MR)

  # extra descriptive stats (relative abundance)
  rel <- sweep(counts, 2, colSums(counts), "/")
  mean_relab_GR <- rowMeans(rel[, s_GR, drop = FALSE])
  mean_relab_MR <- rowMeans(rel[, s_MR, drop = FALSE])

  prev_all <- rowMeans(counts > 0)
  prev_GR  <- rowMeans(counts[, s_GR, drop = FALSE] > 0)
  prev_MR  <- rowMeans(counts[, s_MR, drop = FALSE] > 0)

  all_tbl <- tibble(
    level = lv,
    taxon = rownames(counts),

    delta_CLR = as.numeric(delta_clr),
    log2FC_rel = as.numeric(log2fc_rel),

    mean_relab_GR = as.numeric(mean_relab_GR),
    mean_relab_MR = as.numeric(mean_relab_MR),

    prev_all = as.numeric(prev_all),
    prev_GR  = as.numeric(prev_GR),
    prev_MR  = as.numeric(prev_MR),

    p_value = as.numeric(pvals),
    q_value = as.numeric(qvals),

    direction = case_when(
      delta_CLR > 0 ~ "Higher_in_GR",
      delta_CLR < 0 ~ "Lower_in_GR",
      TRUE ~ "Zero"
    )
  ) %>%
    arrange(q_value, desc(abs(delta_CLR)))

  # outputs
  out_all <- file.path(OUTDIR, paste0("ALL_", lv, "_GR_vs_MR_CLR_Wilcox_BH.csv"))
  write_csv(all_tbl, out_all)

  sig_tbl <- all_tbl %>% filter(!is.na(q_value), q_value < Q_CUT)
  out_sig <- file.path(OUTDIR, paste0("SIGNATURE_q0.1_", lv, "_GR_vs_MR.csv"))
  write_csv(sig_tbl, out_sig)

  summary_rows[[lv]] <- tibble(
    level = lv,
    input_file = obj$file,
    n_tested = nrow(all_tbl),
    n_signature_q0.1 = nrow(sig_tbl),
    q_min = suppressWarnings(min(all_tbl$q_value, na.rm = TRUE))
  )

  message(sprintf("[OK] %s: tested=%d, signature(q<0.1)=%d, q_min=%.4g",
                  lv, nrow(all_tbl), nrow(sig_tbl),
                  suppressWarnings(min(all_tbl$q_value, na.rm = TRUE))))
}

summary_tbl <- bind_rows(summary_rows)
write_csv(summary_tbl, file.path(OUTDIR, "SUMMARY_signature_q0.1_6levels.csv"))

message("\n[DONE] All outputs in folder: ", OUTDIR)
