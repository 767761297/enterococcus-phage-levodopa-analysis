# ============================================================
# PCA (All metabolites): uses all metabolites in raw table
# Ellipse per group (shaded)
# Input (in working directory):
#   定量结果2.xlsx (sheet: 样本浓度)
# Output:
#   PCA_AllMet_ellipse.png / .pdf
# ============================================================

# 1) 安装/加载包
pkgs <- c("readxl","dplyr","stringr","ggplot2","ggrepel")
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# 2) 参数
raw_file  <- "定量结果2.xlsx"
raw_sheet <- "样本浓度"

ELLIPSE_LEVEL <- 0.68     # 0.68≈1 SD；可改为0.95
LABEL_SAMPLES <- FALSE #TURE：显示样本名称

out_png <- "PCA_AllMet_ellipse.png"
out_pdf <- "PCA_AllMet_ellipse.pdf"

# 3) 读取数据
raw <- read_excel(raw_file, sheet = raw_sheet)

# 4) 识别样本列
sample_cols <- names(raw)[str_detect(names(raw), "^[ABCD]\\d+$")]
if (length(sample_cols) == 0) stop("未识别到样本列（A1..D5）。请检查Excel列名。")

# 5) 识别代谢物列（全代谢物 PCA 不依赖 stats 表）
met_col_raw <- intersect(names(raw),
                         c("metabolite","Compound Name (nmol/g)","Compound Name",
                           "代谢物","Metabolite","Name"))[1]
if (is.na(met_col_raw)) met_col_raw <- names(raw)[1]

raw2 <- raw %>%
  mutate(metabolite = as.character(.data[[met_col_raw]])) %>%
  select(metabolite, all_of(sample_cols))

# 6) 0值替代 + log2 + 行z-score
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

# 去掉方差为0或含NA的代谢物（稳健）
keep <- apply(X, 1, function(x) all(is.finite(x)) && sd(x) > 0)
X <- X[keep, , drop = FALSE]

Xz <- t(scale(t(X)))
Xz[!is.finite(Xz)] <- 0

# 7) PCA
pca <- prcomp(t(Xz), center = FALSE, scale. = FALSE)
var_exp <- (pca$sdev^2) / sum(pca$sdev^2)

pc <- as.data.frame(pca$x[, 1:2])
pc$sample <- rownames(pc)
pc$Group <- factor(dplyr::recode(substr(pc$sample, 1, 1),
                                 "A"="GR", "B"="MR", "C"="MR+Phage cocktail", "D"="MR+AFMT"),
                   levels = c("GR","MR","MR+Phage cocktail","MR+AFMT"))

pc1_lab <- paste0("PC1 (", round(var_exp[1]*100, 1), "%)")
pc2_lab <- paste0("PC2 (", round(var_exp[2]*100, 1), "%)")

# 8) 作图（椭圆阴影）
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

# 9) 导出
ggsave(out_png, plot = p, width = 7.2, height = 5.2, dpi = 300)
ggsave(out_pdf, plot = p, width = 7.2, height = 5.2)
cat("Saved:", out_png, "\nSaved:", out_pdf, "\n")