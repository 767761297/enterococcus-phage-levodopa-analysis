## =========================================================
## Paper-style PCoA (Bray-Curtis & Jaccard) from phyloseq RDS
## + betadisper (per-sample + pairwise Wilcoxon, 3 pairs)
## + PERMANOVA (overall + pairwise, 3 pairs)
## 3-group version: GR / MR / MR+Phage cocktail
## Output folder: ./BetaPCoA_outputs_phyloseq_3groups/
## =========================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

rds_file <- "phyloseq_16S.rds"
out_dir  <- "BetaPCoA_outputs_phyloseq_3groups"
if (!file.exists(rds_file)) stop("Missing file: ", rds_file)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ps <- readRDS(rds_file)
if (!inherits(ps, "phyloseq")) stop("Input RDS is not a phyloseq object.")

## ---------- 1) OTU/ASV table ----------
otu <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) otu <- t(otu)
otu <- otu[rowSums(otu) > 0, , drop = FALSE]
if (nrow(otu) < 2 || ncol(otu) < 3) stop("Too few taxa or samples after filtering.")

## ---------- 2) metadata ----------
if (is.null(sample_data(ps, errorIfNULL = FALSE))) {
  stop("No sample_data found in the phyloseq object. Please add metadata before running.")
}

sd <- tryCatch(
  as(sample_data(ps), "data.frame"),
  error = function(e1) {
    tryCatch(
      data.frame(sample_data(ps)),
      error = function(e2) {
        stop(
          "Failed to convert sample_data to data.frame.\n",
          "Try: sd <- as(sample_data(ps), 'data.frame')\n",
          "Or check your phyloseq version. Original errors:\n",
          e1$message, "\n", e2$message
        )
      }
    )
  }
)

sd$Sample <- rownames(sd)

## auto-detect group column
cand <- colnames(sd)
cand_l <- tolower(cand)
preferred <- c("group", "grp", "treatment", "condition", "class", "phenotype")
hit <- cand[match(preferred, cand_l)]
hit <- hit[!is.na(hit)]
if (length(hit) == 0) {
  writeLines(
    paste("Available sample_data columns:\n", paste(cand, collapse = ", ")),
    con = file.path(out_dir, "sampledata_columns.txt")
  )
  stop("Cannot find a group column. See: ", file.path(out_dir, "sampledata_columns.txt"))
}
group_col <- hit[1]
sd$Group_raw <- trimws(as.character(sd[[group_col]]))

## ---------- 3) normalize group labels ----------
normalize_group <- function(x) {
  if (is.na(x) || x == "") return(NA_character_)
  xx <- tolower(trimws(x))
  xx <- gsub("[ _\\-]+", "", xx)

  if (grepl("afmt", xx)) return("MR+AFMT")
  if (grepl("phage", xx) || grepl("cocktail", xx)) return("MR+Phage cocktail")

  if (xx == "gr") return("GR")
  if (xx == "mr") return("MR")

  if (xx %in% c("mrafmt", "mr+afmt")) return("MR+AFMT")
  if (xx %in% c("mrphagecocktail", "mr+phagecocktail", "mrphage", "mr+phage", "mrcocktail", "mr+cocktail")) {
    return("MR+Phage cocktail")
  }

  if (x %in% c("GR", "MR", "MR+AFMT", "MR+Phage cocktail")) return(x)
  return(x)
}

sd$Group <- vapply(sd$Group_raw, normalize_group, character(1))

map_summary <- sd %>% count(Group_raw, Group, sort = TRUE)
write.csv(map_summary, file.path(out_dir, "group_mapping_summary.csv"), row.names = FALSE)

## keep only samples present in OTU table
sd <- sd %>% filter(Sample %in% colnames(otu))
otu <- otu[, sd$Sample, drop = FALSE]

## keep ONLY the 3 expected groups
keep_levels <- c("GR", "MR", "MR+Phage cocktail")
drop_other <- !(sd$Group %in% keep_levels)
if (any(drop_other)) {
  dropped <- sd[drop_other, c("Sample", "Group_raw", "Group"), drop = FALSE]
  write.csv(dropped, file.path(out_dir, "dropped_samples_other_groups.csv"), row.names = FALSE)
  sd <- sd[!drop_other, , drop = FALSE]
  otu <- otu[, sd$Sample, drop = FALSE]
}

drop_na <- is.na(sd$Group)
if (any(drop_na)) {
  dropped <- sd[drop_na, c("Sample", "Group_raw", "Group"), drop = FALSE]
  write.csv(dropped, file.path(out_dir, "dropped_samples_missing_group.csv"), row.names = FALSE)
  sd <- sd[!drop_na, , drop = FALSE]
  otu <- otu[, sd$Sample, drop = FALSE]
}

sd$Group <- factor(sd$Group, levels = keep_levels)
sd$Group <- droplevels(sd$Group)

if (nlevels(sd$Group) < 2) stop("After cleaning, <2 groups remain. Check group_mapping_summary.csv")

## ---------- 4) distances ----------
rel <- sweep(otu, 2, colSums(otu), "/")
rel[is.na(rel)] <- 0

bray <- vegdist(t(rel), method = "bray")
pa   <- (otu > 0) * 1
jacc <- vegdist(t(pa), method = "jaccard", binary = TRUE)

saveRDS(bray, file.path(out_dir, "distance_bray.rds"))
saveRDS(jacc, file.path(out_dir, "distance_jaccard.rds"))

## ---------- 5) PCoA ----------
run_pcoa <- function(d) {
  cs <- cmdscale(d, k = 2, eig = TRUE)
  eig <- cs$eig
  eig_pos <- eig[eig > 0]
  if (length(eig_pos) == 0) {
    var <- c(NA_real_, NA_real_)
  } else {
    var <- eig_pos / sum(eig_pos)
    if (length(var) < 2) var <- c(var, rep(NA_real_, 2 - length(var)))
  }
  coords <- as.data.frame(cs$points)
  colnames(coords) <- c("PCoA1", "PCoA2")
  coords$Sample <- rownames(coords)
  list(coords = coords, var = var[1:2])
}

pcoa_bray <- run_pcoa(bray)
pcoa_jacc <- run_pcoa(jacc)

df_bray <- pcoa_bray$coords %>% left_join(sd %>% select(Sample, Group), by = "Sample")
df_jacc <- pcoa_jacc$coords %>% left_join(sd %>% select(Sample, Group), by = "Sample")

write.csv(df_bray, file.path(out_dir, "PCoA_coords_bray.csv"), row.names = FALSE)
write.csv(df_jacc, file.path(out_dir, "PCoA_coords_jaccard.csv"), row.names = FALSE)

## ---------- 6) plot style ----------
group_colors <- c(
  "GR" = "#7CAE00",
  "MR" = "#F8766D",
  "MR+Phage cocktail" = "#00BFC4"
)

paper_pcoa_plot <- function(df, var2) {
  xlab <- paste0("PCoA1 (", ifelse(is.na(var2[1]), "NA", sprintf("%.1f", 100 * var2[1])), "%)")
  ylab <- paste0("PCoA2 (", ifelse(is.na(var2[2]), "NA", sprintf("%.1f", 100 * var2[2])), "%)")

  df$Group <- factor(df$Group, levels = keep_levels)
  df_ell <- df %>% group_by(Group) %>% filter(n() >= 3) %>% ungroup()

  ggplot(df, aes(PCoA1, PCoA2, color = Group, fill = Group)) +
    stat_ellipse(
      data = df_ell,
      geom = "polygon", type = "t", level = 0.68,
      alpha = 0.18, linewidth = 0
    ) +
    geom_point(size = 3, alpha = 0.95) +
    scale_color_manual(values = group_colors, drop = FALSE) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    labs(x = xlab, y = ylab, color = "Group", fill = "Group") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.6),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.4),
      legend.position = "right",
      axis.title = element_text(size = 10),
      axis.text  = element_text(size = 10)
    )
}

p_bray <- paper_pcoa_plot(df_bray, pcoa_bray$var)
p_jacc <- paper_pcoa_plot(df_jacc, pcoa_jacc$var)

ggsave(file.path(out_dir, "PCoA_BrayCurtis_paperStyle.pdf"), p_bray, width = 6, height = 5)
ggsave(file.path(out_dir, "PCoA_BrayCurtis_paperStyle.png"), p_bray, width = 6, height = 5, dpi = 300)
ggsave(file.path(out_dir, "PCoA_Jaccard_paperStyle.pdf"), p_jacc, width = 6, height = 5)
ggsave(file.path(out_dir, "PCoA_Jaccard_paperStyle.png"), p_jacc, width = 6, height = 5, dpi = 300)

## ---------- 7) betadisper per-sample + pairwise Wilcoxon ----------
pairwise_wilcox <- function(per_sample) {
  lv <- levels(droplevels(per_sample$Group))
  comb <- combn(lv, 2, simplify = FALSE)
  res <- lapply(comb, function(p) {
    d1 <- per_sample$Distance[per_sample$Group == p[1]]
    d2 <- per_sample$Distance[per_sample$Group == p[2]]
    wt <- wilcox.test(d1, d2, exact = FALSE)
    data.frame(
      group1 = p[1], group2 = p[2],
      n1 = length(d1), n2 = length(d2),
      median1 = median(d1, na.rm = TRUE),
      median2 = median(d2, na.rm = TRUE),
      delta_median = median(d1, na.rm = TRUE) - median(d2, na.rm = TRUE),
      stat_W = unname(wt$statistic),
      p_value = wt$p.value
    )
  }) %>% bind_rows()
  res$q_value <- p.adjust(res$p_value, method = "BH")
  res
}

run_betadisper <- function(dist_obj, tag) {
  bd <- betadisper(dist_obj, sd$Group, type = "centroid")

  per_sample <- data.frame(
    Sample   = names(bd$distances),
    Group    = sd$Group[match(names(bd$distances), sd$Sample)],
    Distance = as.numeric(bd$distances)
  )
  per_sample$Group <- droplevels(per_sample$Group)

  write.csv(per_sample, file.path(out_dir, paste0("betadisper_perSample_", tag, ".csv")), row.names = FALSE)
  pw <- pairwise_wilcox(per_sample)
  write.csv(pw, file.path(out_dir, paste0("betadisper_pairwise_Wilcoxon_3pairs_", tag, ".csv")), row.names = FALSE)

  bd_aov <- anova(bd)
  bd_perm <- permutest(bd, permutations = 9999)
  capture.output(bd_aov, bd_perm, file = file.path(out_dir, paste0("betadisper_global_", tag, ".txt")))
}

set.seed(1)
run_betadisper(bray, "bray")
run_betadisper(jacc, "jaccard")

## ---------- 8) PERMANOVA overall + pairwise ----------
run_permanova <- function(dist_obj, tag, permutations = 9999) {
  set.seed(1)
  overall <- adonis2(dist_obj ~ Group, data = sd, permutations = permutations, by = "margin")
  overall_df <- as.data.frame(overall)
  overall_df$term <- rownames(overall_df)
  rownames(overall_df) <- NULL
  write.csv(overall_df, file.path(out_dir, paste0("PERMANOVA_overall_", tag, ".csv")), row.names = FALSE)

  lv <- levels(droplevels(sd$Group))
  comb <- combn(lv, 2, simplify = FALSE)

  pw <- lapply(comb, function(p) {
    keep <- sd$Group %in% p
    sub_sd <- sd[keep, , drop = FALSE]
    sub_dist <- as.dist(as.matrix(dist_obj)[sub_sd$Sample, sub_sd$Sample])

    set.seed(1)
    ad <- adonis2(sub_dist ~ Group, data = sub_sd, permutations = permutations, by = "margin")
    ad_df <- as.data.frame(ad)
    ad_df$term <- rownames(ad_df)
    rownames(ad_df) <- NULL
    g_row <- ad_df %>% filter(term == "Group") %>% slice(1)

    data.frame(
      group1 = p[1], group2 = p[2],
      Df = g_row$Df, R2 = g_row$R2, F = g_row$F,
      p_value = g_row$`Pr(>F)`
    )
  }) %>% bind_rows()

  pw$q_value <- p.adjust(pw$p_value, method = "BH")
  write.csv(pw, file.path(out_dir, paste0("PERMANOVA_pairwise_3pairs_", tag, ".csv")), row.names = FALSE)
}

run_permanova(bray, "bray")
run_permanova(jacc, "jaccard")

message("Done. Outputs in: ", normalizePath(out_dir))
message("Group column used: ", group_col)
message("Groups kept: ", paste(levels(sd$Group), collapse = ", "))