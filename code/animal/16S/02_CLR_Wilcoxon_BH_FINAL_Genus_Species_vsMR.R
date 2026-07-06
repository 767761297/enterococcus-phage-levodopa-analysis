## 02_CLR_Wilcoxon_BH_FINAL_Genus_Species_vsMR.R
## Differential abundance for 16S: CLR + Wilcoxon + BH (FDR)
## - Works on phyloseq_16S.rds
## - Aggregates to Genus + Species
## - 3 planned contrasts vs MR: GR vs MR, MR_Phage vs MR, MR_AFMT vs MR
## - No fixed q-threshold; outputs p_value and q_value for flexible downstream filtering
## - Includes prevalence filtering + NA/constant handling for stability

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(readr)
})

# -------------------- User settings --------------------
PS_RDS    <- "phyloseq_16S.rds"
OUTDIR    <- "16S_CLR_Wilcox_BH_FINAL"
PSEUDO    <- 1           # pseudocount for CLR
PREV_CUT  <- 0.10        # keep taxa present (>0) in >=10% samples
# -------------------------------------------------------

dir.create(OUTDIR, showWarnings = FALSE)

stopifnot(file.exists(PS_RDS))
ps <- readRDS(PS_RDS)

# -------------------- 0) group_safe --------------------
# Expect sample_data(ps)$group exists; otherwise infer from sample names A/B/C/D.
meta <- data.frame(sample_data(ps)) %>% rownames_to_column("sample")

if (!("group" %in% names(meta))) {
  group_of <- function(x){
    if (startsWith(x, "A")) return("GR")
    if (startsWith(x, "B")) return("MR")
    if (startsWith(x, "C")) return("MR+Phage cocktail")
    if (startsWith(x, "D")) return("MR+AFMT")
    return(NA_character_)
  }
  meta$group <- vapply(meta$sample, group_of, character(1))
}

meta$group_safe <- dplyr::case_when(
  meta$group == "MR" ~ "MR",
  meta$group == "GR" ~ "GR",
  str_detect(meta$group, "Phage") ~ "MR_Phage",
  str_detect(meta$group, "AFMT")  ~ "MR_AFMT",
  TRUE ~ NA_character_
)
stopifnot(!any(is.na(meta$group_safe)))

meta$group_safe <- factor(meta$group_safe, levels = c("MR","GR","MR_Phage","MR_AFMT"))

# write back
sd <- data.frame(sample_data(ps)) %>%
  rownames_to_column("sample") %>%
  select(-any_of("group_safe")) %>%
  left_join(meta %>% select(sample, group_safe), by = "sample") %>%
  column_to_rownames("sample")
sample_data(ps) <- sample_data(sd)

# planned contrasts vs MR
CONTRASTS <- list(
  GR_vs_MR       = c("GR", "MR"),
  MR_Phage_vs_MR = c("MR_Phage", "MR"),
  MR_AFMT_vs_MR  = c("MR_AFMT", "MR")
)

# -------------------- helpers --------------------
make_rank_object <- function(ps_in, rank = "Genus") {
  ps_glom <- tax_glom(ps_in, taxrank = rank, NArm = FALSE)
  tax_df <- as.data.frame(tax_table(ps_glom))
  stopifnot(rank %in% colnames(tax_df))

  if (rank == "Species") {
    genus <- if ("Genus" %in% colnames(tax_df)) as.character(tax_df$Genus) else NA_character_
    sp <- as.character(tax_df$Species)
    sp2 <- ifelse(!is.na(sp) & sp != "", sp,
                  ifelse(!is.na(genus) & genus != "", paste0(genus, " sp."), "Unclassified"))
    tax_df$Species <- sp2
    label <- tax_df$Species
  } else {
    tax_df[[rank]] <- ifelse(is.na(tax_df[[rank]]) | tax_df[[rank]] == "", "Unclassified", as.character(tax_df[[rank]]))
    label <- tax_df[[rank]]
  }

  tax_table(ps_glom) <- tax_table(as.matrix(tax_df))
  taxa_names(ps_glom) <- make.unique(label)  # readable unique names
  ps_glom
}

get_counts_taxa_by_sample <- function(ps_in) {
  otu <- as(otu_table(ps_in), "matrix")
  if (!taxa_are_rows(ps_in)) otu <- t(otu)
  otu <- round(otu)
  storage.mode(otu) <- "integer"
  otu
}

prevalence_filter <- function(counts, prev_cut = 0.10) {
  prv <- rowMeans(counts > 0)
  keep <- prv >= prev_cut
  counts[keep, , drop = FALSE]
}

clr_transform <- function(counts, pseudo = 1) {
  x <- counts + pseudo
  logx <- log(x)
  sweep(logx, 2, colMeans(logx), "-")  # CLR per sample
}

# compute log2FC based on mean relative abundance (with tiny pseudo to avoid /0)
log2fc_rel <- function(counts, s1, s0, pseudo_rel = 1e-6) {
  rel <- sweep(counts, 2, colSums(counts), "/")
  m1 <- rowMeans(rel[, s1, drop = FALSE]) + pseudo_rel
  m0 <- rowMeans(rel[, s0, drop = FALSE]) + pseudo_rel
  log2(m1 / m0)
}

safe_wilcox <- function(x1, x0) {
  # if constant across both groups -> p=1
  if (length(unique(c(x1, x0))) <= 1) return(1)
  out <- tryCatch(
    suppressWarnings(wilcox.test(x1, x0, exact = FALSE)$p.value),
    error = function(e) NA_real_
  )
  out
}

run_rank <- function(ps_in, rank_label) {
  ps_r <- make_rank_object(ps_in, rank_label)
  counts <- get_counts_taxa_by_sample(ps_r)
  counts <- prevalence_filter(counts, PREV_CUT)

  # align meta to matrix columns
  meta2 <- data.frame(sample_data(ps_r)) %>%
    rownames_to_column("sample") %>%
    mutate(group_safe = factor(group_safe, levels = c("MR","GR","MR_Phage","MR_AFMT")))

  stopifnot(all(colnames(counts) %in% meta2$sample))
  meta2 <- meta2[match(colnames(counts), meta2$sample), ]
  stopifnot(all(meta2$sample == colnames(counts)))

  clr <- clr_transform(counts, PSEUDO)

  res_list <- list()

  for (nm in names(CONTRASTS)) {
    g1 <- CONTRASTS[[nm]][1]
    g0 <- CONTRASTS[[nm]][2]

    s1 <- meta2$sample[meta2$group_safe == g1]
    s0 <- meta2$sample[meta2$group_safe == g0]

    # effect sizes
    d_clr <- rowMeans(clr[, s1, drop = FALSE]) - rowMeans(clr[, s0, drop = FALSE])
    l2fc  <- log2fc_rel(counts, s1, s0)

    # p-values
    pvals <- apply(clr, 1, function(v) safe_wilcox(v[s1], v[s0]))
    qvals <- p.adjust(pvals, method = "BH")

    df <- tibble(
      rank = rank_label,
      taxon = rownames(counts),
      comparison = nm,
      group1 = g1,
      group0 = g0,
      n_group1 = length(s1),
      n_group0 = length(s0),

      delta_CLR = as.numeric(d_clr),      # mean(CLR_g1) - mean(CLR_g0)
      log2FC_rel = as.numeric(l2fc),      # mean(relative abundance) ratio in log2

      p_value = as.numeric(pvals),
      q_value = as.numeric(qvals)
    ) %>%
      mutate(
        direction = case_when(
          !is.na(q_value) & delta_CLR > 0 ~ "UP_in_group1",
          !is.na(q_value) & delta_CLR < 0 ~ "DOWN_in_group1",
          TRUE ~ NA_character_
        )
      ) %>%
      # drop impossible rows (rare NaNs)
      filter(is.finite(delta_CLR), is.finite(log2FC_rel)) %>%
      arrange(q_value, desc(abs(delta_CLR)))

    res_list[[nm]] <- df
  }

  bind_rows(res_list)
}

# -------------------- run Genus + Species --------------------
cat("[RUN] Genus...\n")
genus_tbl <- run_rank(ps, "Genus")
write_excel_csv(genus_tbl, file.path(OUTDIR, "CLR_WilcoxBH_Genus_ALL.csv"))

cat("[RUN] Species...\n")
species_tbl <- run_rank(ps, "Species")
write_excel_csv(species_tbl, file.path(OUTDIR, "CLR_WilcoxBH_Species_ALL.csv"))

# -------------------- split outputs per comparison --------------------
split_write <- function(tbl, rank_label) {
  for (cmp in unique(tbl$comparison)) {
    fn <- file.path(OUTDIR, paste0("CLR_WilcoxBH_", rank_label, "_", cmp, "_ALL.csv"))
    write_excel_csv(tbl %>% filter(comparison == cmp), fn)
  }
}

split_write(genus_tbl, "Genus")
split_write(species_tbl, "Species")

# -------------------- summary tables (no threshold, just helpful counts) --------------------
summarize_tbl <- function(tbl, rank_label) {
  s <- tbl %>%
    group_by(comparison) %>%
    summarise(
      n_taxa_tested = n(),
      n_na_p = sum(is.na(p_value)),
      n_na_q = sum(is.na(q_value)),
      q_min  = suppressWarnings(min(q_value, na.rm = TRUE)),
      q_p10  = suppressWarnings(quantile(q_value, 0.10, na.rm = TRUE)),
      q_p25  = suppressWarnings(quantile(q_value, 0.25, na.rm = TRUE)),
      q_med  = suppressWarnings(median(q_value, na.rm = TRUE)),
      .groups = "drop"
    )
  write_excel_csv(s, file.path(OUTDIR, paste0("SUMMARY_", rank_label, ".csv")))
}

summarize_tbl(genus_tbl, "Genus")
summarize_tbl(species_tbl, "Species")

# -------------------- session log + parameters --------------------
params <- tibble(
  PS_RDS = PS_RDS,
  PSEUDO = PSEUDO,
  PREV_CUT = PREV_CUT
)
write_excel_csv(params, file.path(OUTDIR, "RUN_parameters.csv"))

sink(file.path(OUTDIR, "SessionInfo.txt"))
sessionInfo()
sink()

cat("[OK] Done. Outputs in folder: ", OUTDIR, "\n")
