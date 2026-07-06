# ============================================================
# 01_diversity_export_tables_FULL6.R
# Export rich diversity tables for local stats/plotting
# - QC read depth per sample
# - Alpha diversity per sample (raw + rarefied), plus:
#     * overall Kruskal-Wallis + BH
#     * pairwise Wilcoxon full 6 comparisons: stats + p + BH-q
# - Beta diversity (Bray/Jaccard/UniFrac/WUniFrac if tree exists):
#     * distance matrix
#     * PCoA coords (per sample) + eigen/variance
#     * PERMANOVA global
#     * PERMANOVA pairwise full 6 comparisons: F/R2 + p + BH-q
#     * betadisper global
#     * betadisper distance-to-centroid per sample
#     * betadisper pairwise full 6 comparisons: stats + p + BH-q
# - Feature tables (counts + rel abundance) and taxonomy (if any)
# ============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

# =========================
# 0) PARAMETERS
# =========================
PS_RDS  <- "phyloseq_16S.rds"
OUTDIR  <- "16S_diversity_exports_FULL6"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# group column priority (auto-detect)
GROUP_COL_PRIORITY <- c("group_safe", "group", "Group", "GROUP", "treatment", "Treatment")

# beta distance methods to compute (UniFrac requires phy_tree)
BETA_METHODS <- c("bray", "jaccard", "unifrac", "wunifrac")

# rarefaction control
DO_RAREFY <- TRUE
RAREFACTION_SEED <- 1
# if you want to force a specific depth instead of min depth, set a number >0
RAREFACTION_DEPTH_OVERRIDE <- 0

# PERMANOVA/betadisper permutations
N_PERM <- 9999

# -------------------------
# Helper: pick group column
# -------------------------
pick_group_col <- function(meta_df) {
  hit <- intersect(GROUP_COL_PRIORITY, names(meta_df))
  if (length(hit) == 0) {
    stop("Cannot find group column in sample_data(ps). Available columns: ",
         paste(names(meta_df), collapse = ", "))
  }
  hit[1]
}

# -------------------------
# Helper: safe CSV write
# -------------------------
safe_write_csv <- function(df, path) {
  write_csv(df, path)
  invisible(TRUE)
}

# =========================
# 1) LOAD DATA + METADATA
# =========================
stopifnot(file.exists(PS_RDS))
ps <- readRDS(PS_RDS)

meta <- data.frame(sample_data(ps)) %>%
  rownames_to_column("sample")

grp_col <- pick_group_col(meta)

meta <- meta %>%
  mutate(
    group = as.factor(.data[[grp_col]]),
    sample = as.character(sample)
  ) %>%
  select(sample, group, everything())

# ensure group levels are stable (A/B/C/D or GR/MR/..)
meta$group <- droplevels(meta$group)

safe_write_csv(meta, file.path(OUTDIR, "Sample_Metadata.csv"))

# =========================
# 2) QC: READ DEPTH PER SAMPLE
# =========================
depth <- sample_sums(ps)
df_depth <- tibble(sample = names(depth), depth = as.numeric(depth)) %>%
  left_join(meta, by = "sample") %>%
  arrange(group, desc(depth))

safe_write_csv(df_depth, file.path(OUTDIR, "QC_ReadDepth_perSample.csv"))

# optional QC barplot (saved to file; does not require Plots tab)
p_depth <- ggplot(df_depth, aes(x = sample, y = depth, fill = group)) +
  geom_col() +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = NULL, y = "Reads per sample")
ggsave(file.path(OUTDIR, "QC_ReadDepth.pdf"), p_depth, width = 9, height = 3.8)

# =========================
# 3) ALPHA DIVERSITY
# =========================
calc_alpha <- function(ps_obj, label = "raw") {
  alpha <- estimate_richness(
    ps_obj,
    measures = c("Observed","Chao1","ACE","Shannon","Simpson","InvSimpson","Fisher")
  ) %>%
    rownames_to_column("sample") %>%
    left_join(meta, by = "sample") %>%
    mutate(
      Evenness_Pielou = ifelse(Observed > 0, Shannon / log(Observed), NA_real_),
      alpha_type = label
    )
  alpha
}

alpha_raw <- calc_alpha(ps, label = "raw")
safe_write_csv(alpha_raw, file.path(OUTDIR, "Alpha_perSample_raw.csv"))

alpha_rarefy <- NULL
ps_rar <- NULL

if (DO_RAREFY) {
  set.seed(RAREFACTION_SEED)
  min_depth <- min(sample_sums(ps))

  if (RAREFACTION_DEPTH_OVERRIDE > 0) {
    if (RAREFACTION_DEPTH_OVERRIDE > min_depth) {
      stop("RAREFACTION_DEPTH_OVERRIDE is higher than min sample depth. ",
           "Override=", RAREFACTION_DEPTH_OVERRIDE, " min_depth=", min_depth)
    }
    rar_depth <- RAREFACTION_DEPTH_OVERRIDE
  } else {
    rar_depth <- min_depth
  }

  ps_rar <- suppressWarnings(
    rarefy_even_depth(
      ps,
      sample.size = rar_depth,
      rngseed = RAREFACTION_SEED,
      replace = FALSE,
      verbose = FALSE
    )
  )

  alpha_rarefy <- calc_alpha(ps_rar, label = paste0("rarefy_", rar_depth))
  safe_write_csv(alpha_rarefy, file.path(OUTDIR, "Alpha_perSample_rarefied.csv"))
  safe_write_csv(tibble(rarefaction_depth = rar_depth),
                 file.path(OUTDIR, "Alpha_rarefaction_depth.csv"))
}

# ---- Alpha overall KW + BH (per alpha_type)
export_alpha_overall_kw <- function(alpha_df, out_csv) {
  metrics <- c("Observed","Chao1","ACE","Shannon","Simpson","InvSimpson","Fisher","Evenness_Pielou")
  out <- lapply(metrics, function(m) {
    dfm <- alpha_df %>% filter(!is.na(.data[[m]]))
    p <- kruskal.test(dfm[[m]] ~ dfm$group)$p.value
    tibble(
      alpha_type = unique(dfm$alpha_type),
      metric = m,
      test = "Kruskal-Wallis",
      p_value = p
    )
  }) %>% bind_rows() %>%
    group_by(alpha_type) %>%
    mutate(q_value = p.adjust(p_value, method = "BH")) %>%
    ungroup()

  safe_write_csv(out, out_csv)
  out
}

# ---- Alpha pairwise full 6 (Wilcoxon): stats + p + BH-q
export_alpha_pairwise_full6 <- function(alpha_df, out_csv) {
  alpha_df <- alpha_df %>% mutate(group = droplevels(as.factor(group)))
  metrics <- c("Observed","Chao1","ACE","Shannon","Simpson","InvSimpson","Fisher","Evenness_Pielou")
  grps <- levels(alpha_df$group)
  pairs <- t(combn(grps, 2))

  out <- lapply(metrics, function(m) {
    dfm <- alpha_df %>% filter(!is.na(.data[[m]]))
    res <- bind_rows(apply(pairs, 1, function(x) {
      g1 <- x[1]; g2 <- x[2]
      v1 <- dfm %>% filter(group == g1) %>% pull(.data[[m]])
      v2 <- dfm %>% filter(group == g2) %>% pull(.data[[m]])
      wt <- suppressWarnings(wilcox.test(v1, v2, exact = FALSE))

      tibble(
        alpha_type = unique(dfm$alpha_type),
        metric = m,
        group1 = g1, group2 = g2,
        n1 = length(v1), n2 = length(v2),
        median1 = median(v1, na.rm = TRUE),
        median2 = median(v2, na.rm = TRUE),
        delta_median = median1 - median2,
        stat_W = unname(wt$statistic),
        p_value = wt$p.value
      )
    })) %>%
      mutate(q_value = p.adjust(p_value, method = "BH")) # BH across 6 comparisons within metric
    res
  }) %>% bind_rows()

  safe_write_csv(out, out_csv)
  out
}

export_alpha_overall_kw(alpha_raw, file.path(OUTDIR, "raw_Alpha_overall_KW_BH.csv"))
export_alpha_pairwise_full6(alpha_raw, file.path(OUTDIR, "raw_Alpha_pairwise_Wilcoxon_full6_stats_p_q.csv"))

if (!is.null(alpha_rarefy)) {
  export_alpha_overall_kw(alpha_rarefy, file.path(OUTDIR, "rarefied_Alpha_overall_KW_BH.csv"))
  export_alpha_pairwise_full6(alpha_rarefy, file.path(OUTDIR, "rarefied_Alpha_pairwise_Wilcoxon_full6_stats_p_q.csv"))
}

# Also export a combined alpha table (raw + rarefied) for convenience
alpha_all <- bind_rows(alpha_raw, alpha_rarefy)
if (nrow(alpha_all) > 0) {
  safe_write_csv(alpha_all, file.path(OUTDIR, "Alpha_perSample_raw_plus_rarefied.csv"))
}

# =========================
# 4) BETA DIVERSITY EXPORTS
# =========================
# Common practice: compute Bray/Jaccard on relative abundance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
has_tree <- !is.null(phy_tree(ps, errorIfNULL = FALSE))

# PERMANOVA pairwise full 6: F/R2 + p + BH-q (n1/n2 included)
export_permanova_pairwise_full6 <- function(dist_obj, meta_df, method, out_csv, n_perm = 9999) {
  meta_df <- meta_df %>% mutate(group = droplevels(as.factor(group)))
  grps <- levels(meta_df$group)
  pairs <- t(combn(grps, 2))
  dist_mat <- as.matrix(dist_obj)

  pw <- apply(pairs, 1, function(x) {
    g1 <- x[1]; g2 <- x[2]
    keep <- meta_df$group %in% c(g1, g2)
    samp <- meta_df$sample[keep]
    dist_sub <- as.dist(dist_mat[samp, samp])

    meta_sub <- meta_df %>%
      filter(sample %in% samp) %>%
      select(sample, group) %>%
      column_to_rownames("sample")

    set.seed(1)
    ad2 <- vegan::adonis2(dist_sub ~ group, data = meta_sub, permutations = n_perm)
    ad2_df <- as.data.frame(ad2) %>% rownames_to_column("term")
    g_row <- ad2_df[ad2_df$term == "group", , drop = FALSE]

    tibble(
      method = method,
      group1 = g1, group2 = g2,
      n1 = sum(meta_df$group == g1),
      n2 = sum(meta_df$group == g2),
      F = g_row$F,
      R2 = g_row$R2,
      p_value = g_row$`Pr(>F)`
    )
  }) %>% bind_rows() %>%
    mutate(q_value = p.adjust(p_value, method = "BH")) # BH across 6 comparisons within method

  safe_write_csv(pw, out_csv)
  pw
}

# betadisper pairwise full 6 on distance-to-centroid: stats + p + BH-q
export_betadisper_pairwise_full6 <- function(bd_dist_df, out_csv) {
  bd_dist_df <- bd_dist_df %>% mutate(group = droplevels(as.factor(group)))
  grps <- levels(bd_dist_df$group)
  pairs <- t(combn(grps, 2))

  out <- bind_rows(apply(pairs, 1, function(x) {
    g1 <- x[1]; g2 <- x[2]
    v1 <- bd_dist_df %>% filter(group == g1) %>% pull(distance_to_centroid)
    v2 <- bd_dist_df %>% filter(group == g2) %>% pull(distance_to_centroid)
    wt <- suppressWarnings(wilcox.test(v1, v2, exact = FALSE))

    tibble(
      method = unique(bd_dist_df$method),
      group1 = g1, group2 = g2,
      n1 = length(v1), n2 = length(v2),
      median1 = median(v1, na.rm = TRUE),
      median2 = median(v2, na.rm = TRUE),
      delta_median = median1 - median2,
      stat_W = unname(wt$statistic),
      p_value = wt$p.value
    )
  })) %>%
    mutate(q_value = p.adjust(p_value, method = "BH")) # BH across 6 comparisons within method

  safe_write_csv(out, out_csv)
  out
}

do_one_beta <- function(method) {
  if (method %in% c("unifrac","wunifrac") && !has_tree) {
    message("[SKIP] ", method, ": phy_tree(ps) not found.")
    return(invisible(NULL))
  }

  message("[BETA] method = ", method)

  # 1) distance
  dist_obj <- phyloseq::distance(ps_rel, method = method)
  dist_mat <- as.matrix(dist_obj)

  safe_write_csv(as.data.frame(dist_mat) %>% rownames_to_column("sample"),
                 file.path(OUTDIR, paste0("Beta_distance_", method, "_matrix.csv")))

  # 2) PCoA coords via cmdscale
  k <- 5
  pcoa <- cmdscale(dist_obj, k = k, eig = TRUE, add = TRUE)
  eig <- pcoa$eig

  # percentage explained: only positive eigenvalues contribute
  pos_eig <- eig[eig > 0]
  var_pct <- rep(NA_real_, length(eig))
  if (length(pos_eig) > 0) {
    var_pct[eig > 0] <- eig[eig > 0] / sum(pos_eig) * 100
  }

  pcoa_df <- as.data.frame(pcoa$points) %>%
    rownames_to_column("sample") %>%
    rename_with(~ paste0("PCoA", seq_along(.x)), -sample) %>%
    left_join(meta, by = "sample") %>%
    mutate(method = method)

  safe_write_csv(pcoa_df, file.path(OUTDIR, paste0("PCoA_coords_", method, ".csv")))

  eig_df <- tibble(
    method = method,
    axis = paste0("PCoA", seq_along(eig)),
    eigenvalue = as.numeric(eig),
    var_explained_pct = as.numeric(var_pct)
  )
  safe_write_csv(eig_df, file.path(OUTDIR, paste0("PCoA_eigen_", method, ".csv")))

  # 3) PERMANOVA global
  meta_vegan <- meta %>% select(sample, group) %>% column_to_rownames("sample")
  set.seed(1)
  ad <- vegan::adonis2(dist_obj ~ group, data = meta_vegan, permutations = N_PERM)
  ad_df <- as.data.frame(ad) %>% rownames_to_column("term") %>% mutate(method = method)
  safe_write_csv(ad_df, file.path(OUTDIR, paste0("PERMANOVA_global_", method, ".csv")))

  # 4) PERMANOVA pairwise full 6
  export_permanova_pairwise_full6(
    dist_obj = dist_obj,
    meta_df  = meta,
    method   = method,
    out_csv  = file.path(OUTDIR, paste0("PERMANOVA_pairwise_", method, "_full6_stats_p_q.csv")),
    n_perm   = N_PERM
  )

  # 5) betadisper global + per-sample distances
  bd <- vegan::betadisper(dist_obj, meta_vegan$group)
  bd_aov <- anova(bd)
  bd_perm <- permutest(bd, permutations = N_PERM)

  bd_global <- tibble(
    method = method,
    anova_F = bd_aov$`F value`[1],
    anova_p = bd_aov$`Pr(>F)`[1],
    perm_F  = bd_perm$tab$F[1],
    perm_p  = bd_perm$tab$`Pr(>F)`[1]
  )
  safe_write_csv(bd_global, file.path(OUTDIR, paste0("betadisper_global_", method, ".csv")))

  bd_dist <- tibble(
    sample = names(bd$distances),
    distance_to_centroid = as.numeric(bd$distances),
    method = method
  ) %>%
    left_join(meta, by = "sample") %>%
    arrange(group, distance_to_centroid)

  safe_write_csv(bd_dist, file.path(OUTDIR, paste0("betadisper_distance_to_centroid_", method, "_perSample.csv")))

  # 6) betadisper pairwise full 6 (Wilcoxon on distance-to-centroid)
  export_betadisper_pairwise_full6(
    bd_dist,
    file.path(OUTDIR, paste0("betadisper_pairwise_Wilcoxon_full6_stats_p_q_", method, ".csv"))
  )

  invisible(NULL)
}

for (m in BETA_METHODS) do_one_beta(m)

# =========================
# 5) FEATURE TABLE EXPORTS (counts + rel abundance)
# =========================
otu_counts <- as.data.frame(otu_table(ps))
if (!taxa_are_rows(ps)) otu_counts <- t(otu_counts)
otu_counts <- as.data.frame(otu_counts) %>% rownames_to_column("feature")
safe_write_csv(otu_counts, file.path(OUTDIR, "FeatureTable_counts.csv"))

otu_rel <- as.data.frame(otu_table(ps_rel))
if (!taxa_are_rows(ps_rel)) otu_rel <- t(otu_rel)
otu_rel <- as.data.frame(otu_rel) %>% rownames_to_column("feature")
safe_write_csv(otu_rel, file.path(OUTDIR, "FeatureTable_relAbundance.csv"))

tax <- tryCatch(as.data.frame(tax_table(ps)), error = function(e) NULL)
if (!is.null(tax)) {
  tax <- tax %>% rownames_to_column("feature")
  safe_write_csv(tax, file.path(OUTDIR, "Taxonomy_table.csv"))
}

# =========================
# 6) SESSION INFO (for reproducibility)
# =========================
sink(file.path(OUTDIR, "sessionInfo.txt"))
print(sessionInfo())
sink()

message("[DONE] All exports saved to: ", OUTDIR)
