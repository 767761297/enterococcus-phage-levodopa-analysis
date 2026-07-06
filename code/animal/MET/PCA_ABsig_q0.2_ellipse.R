# ============================================================
# PCA (AB signature): fdr_A_vs_B < 0.2 (fallback to Top 20)
# Ellipse per group (shaded)
# Inputs (in working directory):
#   1) 定量结果2.xlsx (sheet: 样本浓度)
#   2) neurotransmitter_metabolomics_stats_q0_2.xlsx (sheet: all_results)
# Outputs:
#   PCA_ABsig_q0.2_ellipse.png / .pdf
# ============================================================

# 1) 安装/加载包
pkgs <- c("readxl","dplyr","stringr","ggplot2","ggrepel")
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# 2) 参数
raw_file    <- "定量结果2.xlsx"
raw_sheet   <- "样本浓度"
stats_file  <- "neurotransmitter_metabolomics_stats_q0_2.xlsx"
stats_sheet <- "all_results"

FDR_SIG <- 0.2
TOP_N_FALLBACK <- 20
ELLIPSE_LEVEL <- 0.68     # 0.68≈1 SD；可改为0.95
LABEL_SAMPLES <- TRUE #FALSE :不显示样本名称

out_png <- "PCA_ABsig_q0.2_ellipse.png"
out_pdf <- "PCA_ABsig_q0.2_ellipse.pdf"

# 3) 读取数据
raw <- read_excel(raw_file, sheet = raw_sheet)
res <- read_excel(stats_file, sheet = stats_sheet)

# 4) 识别样本列
sample_cols <- names(raw)[str_detect(names(raw), "^[ABCD]\\d+$")]
if (length(sample_cols) == 0) stop("未识别到样本列（A1..D5）。请检查Excel列名。")

# 5) 识别原始表代谢物列（应为英文以匹配统计表 metabolite）
met_col_raw <- intersect(names(raw),
                         c("metabolite","Compound Name (nmol/g)","Compound Name",
                           "代谢物","Metabolite","Name"))[1]
if (is.na(met_col_raw)) met_col_raw <- names(raw)[1]

raw2 <- raw %>%
  mutate(metabolite = as.character(.data[[met_col_raw]])) %>%
  select(metabolite, all_of(sample_cols))

# 6) 选择 AB 指纹代谢物（q<0.2），过少则退回 Top 20
sig_mets <- res %>%
  filter(!is.na(fdr_A_vs_B), fdr_A_vs_B < FDR_SIG) %>%
  arrange(fdr_A_vs_B) %>%
  pull(metabolite) %>%
  unique()

if (length(sig_mets) < 2) {
  sig_mets <- res %>%
    filter(!is.na(fdr_A_vs_B)) %>%
    arrange(fdr_A_vs_B) %>%
    slice(1:TOP_N_FALLBACK) %>%
    pull(metabolite) %>%
    unique()
  cat("AB signature <2, fallback to Top", TOP_N_FALLBACK, "by fdr_A_vs_B.\n")
} else {
  cat("Using AB signature (q<", FDR_SIG, "), n=", length(sig_mets), "\n", sep = "")
}

# 匹配检查
n_match <- sum(sig_mets %in% raw2$metabolite)
if (n_match < length(sig_mets)) {
  stop(paste0("原始表代谢物名称与统计表英文 metabolite 匹配不足：",
              n_match, "/", length(sig_mets),
              "。请确保原始表代谢物列为英文且与统计表一致。"))
}

# 7) 0值替代（每代谢物 min_pos/2）+ log2 + 行z-score
impute_row_halfmin <- function(v){
  v <- as.numeric(v)
  min_pos <- suppressWarnings(min(v[v > 0], na.rm = TRUE))
  if (is.infinite(min_pos) || is.na(min_pos)) return(rep(NA_real_, length(v)))
  v[v == 0] <- min_pos / 2
  v
}

df <- raw2 %>% filter(metabolite %in% sig_mets) %>%
  rowwise() %>%
  mutate(across(all_of(sample_cols), ~ impute_row_halfmin(.x))) %>%
  ungroup() %>%
  mutate(across(all_of(sample_cols), log2))

mat <- df %>% select(metabolite, all_of(sample_cols)) %>% as.data.frame()
rownames(mat) <- make.unique(mat$metabolite)
mat$metabolite <- NULL
X <- as.matrix(mat)

# 去掉方差为0或含NA的代谢物（稳健）
keep <- apply(X, 1, function(x) all(is.finite(x)) && sd(x) > 0)
X <- X[keep, , drop = FALSE]

Xz <- t(scale(t(X)))
Xz[!is.finite(Xz)] <- 0

# 8) PCA（样本为观测）
pca <- prcomp(t(Xz), center = FALSE, scale. = FALSE)
var_exp <- (pca$sdev^2) / sum(pca$sdev^2)

pc <- as.data.frame(pca$x[, 1:2])
pc$sample <- rownames(pc)
pc$Group <- factor(dplyr::recode(substr(pc$sample, 1, 1),
                                "A"="GR", "B"="MR", "C"="MR+Phage cocktail", "D"="MR+AFMT"),
                   levels = c("GR","MR","MR+Phage cocktail","MR+AFMT"))

pc1_lab <- paste0("PC1 (", round(var_exp[1]*100, 1), "%)")
pc2_lab <- paste0("PC2 (", round(var_exp[2]*100, 1), "%)")

# 9) 作图（椭圆阴影）
p <- ggplot(pc, aes(PC1, PC2)) +
  stat_ellipse(aes(fill = Group), geom = "polygon",
               type = "t", level = ELLIPSE_LEVEL, alpha = 0.18, color = NA) +
  stat_ellipse(aes(color = Group), type = "t",
               level = ELLIPSE_LEVEL, linewidth = 0.6) +
  geom_point(aes(color = Group), size = 3) +
  labs(x = pc1_lab, y = pc2_lab) +
  theme_classic() +
  theme(legend.title = element_blank())

if (LABEL_SAMPLES) {
  p <- p + ggrepel::geom_text_repel(aes(label = sample),
                                   size = 3, max.overlaps = 100)
}

# 10) 导出
ggsave(out_png, plot = p, width = 7.2, height = 5.2, dpi = 300)
ggsave(out_pdf, plot = p, width = 7.2, height = 5.2)
cat("Saved:", out_png, "\nSaved:", out_pdf, "\n")
