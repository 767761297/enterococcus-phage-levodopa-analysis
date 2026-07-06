## 01_diversity_QC_alpha_beta.R
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(phyloseq)
  library(ggplot2)
  library(vegan)
})

PS_RDS <- "phyloseq_16S.rds"
OUTDIR <- "16S_figs_diversity"
dir.create(OUTDIR, showWarnings = FALSE)

ps <- readRDS(PS_RDS)
meta <- data.frame(sample_data(ps)) %>% rownames_to_column("sample")

# ---------- (1) 测序深度 ----------
depth <- sample_sums(ps)
df_depth <- tibble(sample = names(depth), depth = as.numeric(depth)) %>%
  left_join(meta, by = "sample")

p_depth <- ggplot(df_depth, aes(x = sample, y = depth, fill = group)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = NULL, y = "Reads per sample")

ggsave(file.path(OUTDIR, "QC_read_depth.pdf"), p_depth, width = 8, height = 3.5)
write.csv(df_depth, file.path(OUTDIR, "QC_read_depth.csv"), row.names = FALSE)

# ---------- (2) Alpha diversity ----------
alpha <- estimate_richness(ps, measures = c("Observed","Shannon","Simpson")) %>%
  rownames_to_column("sample") %>%
  left_join(meta, by="sample")

plot_alpha <- function(metric){
  ggplot(alpha, aes(x = group, y = .data[[metric]], fill = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.65) +
    geom_jitter(width = 0.12, size = 2) +
    theme_bw() +
    labs(x=NULL, y=metric)
}

for (m in c("Observed","Shannon","Simpson")) {
  g <- plot_alpha(m)
  ggsave(file.path(OUTDIR, paste0("Alpha_", m, ".pdf")), g, width = 6.5, height = 4.2)
}

# 统计：Kruskal-Wallis（整体）
kw <- lapply(c("Observed","Shannon","Simpson"), function(m){
  data.frame(
    metric = m,
    p_value = kruskal.test(alpha[[m]] ~ alpha$group)$p.value
  )
}) %>% bind_rows()
write.csv(kw, file.path(OUTDIR, "Alpha_KruskalWallis.csv"), row.names = FALSE)

# ---------- (3) Beta diversity：Bray-Curtis PCoA + 椭圆 ----------
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
dist_bc <- phyloseq::distance(ps_rel, method = "bray")

ord <- ordinate(ps_rel, method = "PCoA", distance = dist_bc)
ord_df <- as.data.frame(ord$vectors[,1:2])
ord_df$sample <- rownames(ord_df)
ord_df <- ord_df %>% left_join(meta, by="sample")

var_expl <- ord$values$Relative_eig[1:2] * 100

p_pcoa <- ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = group, fill = group)) +
  geom_point(size = 3) +
  stat_ellipse(type = "t", level = 0.68, alpha = 0.18, geom = "polygon", color = NA) +
  theme_bw() +
  labs(x = paste0("PCoA1 (", sprintf("%.1f", var_expl[1]), "%)"),
       y = paste0("PCoA2 (", sprintf("%.1f", var_expl[2]), "%)"))

ggsave(file.path(OUTDIR, "PCoA_BrayCurtis_ellipse.pdf"), p_pcoa, width = 6, height = 5)

# ---------- (4) PERMANOVA + betadisper ----------
df_meta_for_vegan <- meta %>% column_to_rownames("sample")
set.seed(1)
adon <- vegan::adonis2(as.matrix(dist_bc) ~ group, data = df_meta_for_vegan, permutations = 9999)
capture.output(adon, file = file.path(OUTDIR, "PERMANOVA_BrayCurtis.txt"))

bd <- betadisper(dist_bc, df_meta_for_vegan$group)
bd_aov <- anova(bd)
bd_perm <- permutest(bd, permutations = 9999)
capture.output(bd_aov, bd_perm, file = file.path(OUTDIR, "betadisper_BrayCurtis.txt"))

cat("[OK] Diversity outputs saved in:", OUTDIR, "\n")