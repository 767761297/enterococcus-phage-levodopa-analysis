# ============================================================
# All-metabolites feature space
# PERMANOVA (adonis2) + betadisper
# Input (working directory):
#   - 定量结果2.xlsx  (sheet: 样本浓度)
# ============================================================

set.seed(1)

pkgs <- c("readxl","dplyr","stringr","ggplot2","vegan")
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- parameters
raw_file  <- "定量结果2.xlsx"
raw_sheet <- "样本浓度"
N_PERM <- 9999

tag <- "AllMet_featureSpace"

# ---- read
raw <- read_excel(raw_file, sheet = raw_sheet)

sample_cols <- names(raw)[str_detect(names(raw), "^[ABCD]\\d+$")]
if (length(sample_cols) == 0) stop("未识别到样本列（A1..D5）。请检查Excel列名。")

met_col_raw <- intersect(
  names(raw),
  c("metabolite","Compound Name (nmol/g)","Compound Name","代谢物","Metabolite","Name")
)[1]
if (is.na(met_col_raw)) met_col_raw <- names(raw)[1]

raw2 <- raw %>%
  mutate(metabolite = as.character(.data[[met_col_raw]])) %>%
  select(metabolite, all_of(sample_cols))

group_of_sample <- function(s){
  dplyr::recode(substr(s, 1, 1),
                "A"="GR", "B"="MR", "C"="MR+Phage cocktail", "D"="MR+AFMT")
}

# ---- preprocess: 0->min_pos/2, log2, row z-score
impute_row_halfmin <- function(v){
  v <- as.numeric(v)
  min_pos <- suppressWarnings(min(v[v > 0], na.rm = TRUE))
  if (is.infinite(min_pos) || is.na(min_pos)) return(rep(NA_real_, length(v)))
  v[v == 0] <- min_pos / 2
  v
}

df <- raw2 %>%
  rowwise() %>%
  mutate(across(all_of(sample_cols), ~ impute_row_halfmin(.x))) %>%
  ungroup() %>%
  mutate(across(all_of(sample_cols), log2))

mat <- df %>% select(metabolite, all_of(sample_cols)) %>% as.data.frame()
rownames(mat) <- make.unique(mat$metabolite)
mat$metabolite <- NULL
X <- as.matrix(mat)

keep <- apply(X, 1, function(x) all(is.finite(x)) && sd(x) > 0)
X <- X[keep, , drop = FALSE]

Xz <- t(scale(t(X)))
Xz[!is.finite(Xz)] <- 0

# ---- build sample metadata + distance
Y <- t(Xz)
meta <- data.frame(
  sample = rownames(Y),
  Group  = factor(vapply(rownames(Y), group_of_sample, character(1)),
                  levels = c("GR","MR","MR+Phage cocktail","MR+AFMT")),
  stringsAsFactors = FALSE
)

d <- dist(Y, method = "euclidean")

# ---- PERMANOVA (global)
ad <- adonis2(d ~ Group, data = meta, permutations = N_PERM)
ad_df <- as.data.frame(ad)
ad_df$term <- rownames(ad_df)
rownames(ad_df) <- NULL
ad_df <- ad_df %>% select(term, everything())
write.csv(ad_df, paste0("PERMANOVA_global_", tag, ".csv"), row.names = FALSE)

# ---- PERMANOVA (pairwise)
lv <- levels(meta$Group)
pairs <- combn(lv, 2, simplify = FALSE)

pair_res <- lapply(pairs, function(p){
  g1 <- p[1]; g2 <- p[2]
  idx <- meta$Group %in% c(g1, g2)
  meta2 <- meta[idx, , drop = FALSE]
  Y2 <- Y[idx, , drop = FALSE]
  d2 <- dist(Y2, method = "euclidean")
  ad2 <- adonis2(d2 ~ Group, data = meta2, permutations = N_PERM)
  a <- as.data.frame(ad2)
  data.frame(group1=g1, group2=g2, R2=a$R2[1], F=a$F[1], p_value=a$`Pr(>F)`[1])
}) %>% bind_rows() %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

write.csv(pair_res, paste0("PERMANOVA_pairwise_", tag, ".csv"), row.names = FALSE)

# ---- betadisper
bd <- betadisper(d, meta$Group)
bd_aov <- anova(bd)
bd_perm <- permutest(bd, permutations = N_PERM)

capture.output(
  list("ANOVA(betadisper):", bd_aov, "", "permutest(betadisper):", bd_perm, "", "TukeyHSD(betadisper):", TukeyHSD(bd)),
  file = paste0("betadisper_", tag, ".txt")
)

disp_df <- data.frame(sample=meta$sample, Group=meta$Group, dist_to_centroid=bd$distances)
write.csv(disp_df, paste0("betadisper_distances_", tag, ".csv"), row.names = FALSE)

p_disp <- ggplot(disp_df, aes(x=Group, y=dist_to_centroid)) +
  geom_boxplot(width=0.55, outlier.shape=NA, linewidth=0.7) +
  geom_jitter(width=0.12, size=2.2, alpha=0.85) +
  theme_classic(base_size=14) +
  labs(x=NULL, y="Distance to group centroid (dispersion)")

ggsave(paste0("betadisper_boxplot_", tag, ".png"), p_disp, width=7.6, height=5.0, dpi=300)
ggsave(paste0("betadisper_boxplot_", tag, ".pdf"), p_disp, width=7.6, height=5.0)

# ---- PCoA plot
pcoa <- cmdscale(d, k=2, eig=TRUE)
pts <- as.data.frame(pcoa$points)
colnames(pts) <- c("Axis1","Axis2")
pts$sample <- rownames(pts)
pts <- pts %>% left_join(meta, by="sample")

var1 <- round(100 * pcoa$eig[1] / sum(pcoa$eig), 1)
var2 <- round(100 * pcoa$eig[2] / sum(pcoa$eig), 1)

p_pcoa <- ggplot(pts, aes(Axis1, Axis2, color=Group, fill=Group)) +
  stat_ellipse(geom="polygon", type="t", level=0.68, alpha=0.18, color=NA) +
  stat_ellipse(type="t", level=0.68, linewidth=0.6) +
  geom_point(size=3) +
  theme_classic(base_size=14) +
  theme(legend.title=element_blank()) +
  labs(x=paste0("PCoA1 (", var1, "%)"), y=paste0("PCoA2 (", var2, "%)"))

ggsave(paste0("PCoA_", tag, ".png"), p_pcoa, width=7.2, height=5.2, dpi=300)
ggsave(paste0("PCoA_", tag, ".pdf"), p_pcoa, width=7.2, height=5.2)

cat("DONE:", tag, "\nOutputs saved in:\n", getwd(), "\n")
