## =========================================================
## STEP_reduce_corr_heatmap_from_CORR_EDGES.R
## 用两个输入表（CORR + EDGES）把相关热图“瘦身”
## 输出：精简后的 taxa/edge 表 + 精简热图 PDF/PNG
## =========================================================

## ---------- 0) 参数区：你只需要改这里 ----------
CORR_FILE  <- "CORR_species_16Sfp_q0.10__Metfp_q0.20_Spearman.csv"
EDGES_FILE <- "EDGES_species_q0.10_rho0.50.csv"

OUTDIR <- "Corr_heatmap_reduced_species"

## 精简策略（三选一）：
## 1) "topK_per_met"  : 每个代谢物保留 TopK 条关联（推荐；保证每列都有代表）
## 2) "topN_edges"    : 全局保留 TopN 条关联（最狠，可能丢一些代谢物列）
## 3) "topN_taxa_deg" : 按 taxa 的“连接数/度”选 TopN taxa，再画它们与全部代谢物的热图
REDUCE_METHOD <- "topN_edges"

topK_per_met <- 2      # 方法1：每个代谢物保留 K 条（建议 1~3）
topN_edges   <- 10     # 方法2：全局保留 N 条
topN_taxa_deg    <- 20     # 方法3：保留连接数最高的 N 个 taxa

## 强制保留的菌（例如你想文章里一定出现 E. faecalis）
#FORCE_INCLUDE_TAXA <- c("Enterococcus_faecalis")

## 画图细节
PLOT_TITLE <- NULL     # 你说要更干净就不要标题
FONT_SIZE_X <- 10
FONT_SIZE_Y <- 8

## ---------- 1) 依赖包 ----------
req_pkgs <- c("readr","dplyr","tidyr","ggplot2","stringr","forcats")
for (p in req_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
library(readr); library(dplyr); library(tidyr); library(ggplot2)
library(stringr); library(forcats)

## ---------- 2) 读入数据 ----------
stopifnot(file.exists(CORR_FILE), file.exists(EDGES_FILE))

corr  <- readr::read_csv(CORR_FILE,  show_col_types = FALSE)
edges <- readr::read_csv(EDGES_FILE, show_col_types = FALSE)

## 兼容列名（你的文件目前是这些列；这里做个鲁棒防护）
need_corr  <- c("taxon","metabolite","rho","p_value","q_value","n")
need_edges <- c("taxon","metabolite","rho","p_value","q_value","n")

miss_corr  <- setdiff(need_corr,  names(corr))
miss_edges <- setdiff(need_edges, names(edges))
if (length(miss_corr) > 0)  stop("CORR_FILE 缺列: ",  paste(miss_corr, collapse = ", "))
if (length(miss_edges) > 0) stop("EDGES_FILE 缺列: ", paste(miss_edges, collapse = ", "))

corr <- corr %>%
  mutate(
    taxon      = as.character(taxon),
    metabolite = as.character(metabolite),
    rho        = as.numeric(rho),
    p_value    = as.numeric(p_value),
    q_value    = as.numeric(q_value),
    n          = as.integer(n),
    abs_rho    = abs(rho)
  )

edges <- edges %>%
  mutate(
    taxon      = as.character(taxon),
    metabolite = as.character(metabolite),
    rho        = as.numeric(rho),
    p_value    = as.numeric(p_value),
    q_value    = as.numeric(q_value),
    n          = as.integer(n),
    abs_rho    = abs(rho)
  )

## 代谢物集合（通常是 13 个代谢物指纹）
met_set <- sort(unique(edges$metabolite))

## ---------- 3) 精简逻辑：产出 edges_keep + taxa_keep ----------
edges_keep <- NULL
taxa_keep  <- NULL

if (REDUCE_METHOD == "topK_per_met") {
  edges_keep <- edges %>%
    group_by(metabolite) %>%
    arrange(q_value, desc(abs_rho), desc(n)) %>%
    slice_head(n = TOPK_PER_MET) %>%
    ungroup()

} else if (REDUCE_METHOD == "topN_edges") {
  edges_keep <- edges %>%
    arrange(q_value, desc(abs_rho), desc(n)) %>%
    slice_head(n = TOPN_EDGES)

} else if (REDUCE_METHOD == "topN_taxa_deg") {
  taxa_rank <- edges %>%
    count(taxon, name = "degree") %>%
    arrange(desc(degree), taxon)

  taxa_keep <- taxa_rank %>%
    slice_head(n = TOPN_TAXA) %>%
    pull(taxon)

  edges_keep <- edges %>%
    filter(taxon %in% taxa_keep)

} else {
  stop("Unknown REDUCE_METHOD: ", REDUCE_METHOD)
}

## 强制加入指定 taxa（若它在 edges 里出现）
if (length(FORCE_INCLUDE_TAXA) > 0) {
  edges_force <- edges %>% filter(taxon %in% FORCE_INCLUDE_TAXA)
  edges_keep <- bind_rows(edges_keep, edges_force) %>%
    distinct(taxon, metabolite, .keep_all = TRUE)
}

taxa_keep <- sort(unique(edges_keep$taxon))

## ---------- 4) 输出“精简表” ----------
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

## 4.1 精简 edges
edges_keep_out <- edges_keep %>%
  arrange(q_value, desc(abs_rho), taxon, metabolite)

write_csv(edges_keep_out, file.path(OUTDIR, "EDGES_reduced.csv"))

## 4.2 taxa 汇总（连接数/强度/最显著 q）
taxa_summary <- edges %>%
  mutate(in_reduced = taxon %in% taxa_keep) %>%
  group_by(taxon) %>%
  summarise(
    degree_all   = n(),
    degree_in    = sum(in_reduced),
    min_q_all    = min(q_value, na.rm = TRUE),
    max_abs_rho  = max(abs_rho, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(degree_in), min_q_all, desc(max_abs_rho))

write_csv(taxa_summary %>% filter(taxon %in% taxa_keep),
          file.path(OUTDIR, "TAXA_reduced_summary.csv"))

## ---------- 5) 用 CORR 表为“精简 taxa × 全部代谢物集合”补齐 rho，画热图 ----------
plot_df <- corr %>%
  filter(taxon %in% taxa_keep, metabolite %in% met_set) %>%
  select(taxon, metabolite, rho, q_value) %>%
  distinct()

## 如果某些组合在 CORR 中缺失（理论上不该），补 NA
plot_df <- plot_df %>%
  complete(taxon, metabolite, fill = list(rho = NA_real_, q_value = NA_real_))

## 排序：先按 taxa 在 edges_keep 的“连接数”降序；代谢物按字母（你也可改成按生物学顺序）
taxa_order <- edges_keep_out %>% count(taxon, name = "deg") %>% arrange(desc(deg), taxon) %>% pull(taxon)
met_order  <- met_set

plot_df <- plot_df %>%
  mutate(
    taxon      = factor(taxon, levels = taxa_order),
    metabolite = factor(metabolite, levels = met_order)
  )

## 画图（不在格子上写菌名/代谢物名以外的内容；重点是缩短行数）
p <- ggplot(plot_df, aes(x = metabolite, y = taxon, fill = rho)) +
  geom_tile(color = NA) +
  scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#c53030", midpoint = 0,
                       na.value = "grey95", name = "Spearman rho") +
  labs(title = PLOT_TITLE, x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = FONT_SIZE_X),
    axis.text.y = element_text(size = FONT_SIZE_Y),
    plot.title  = element_text(hjust = 0.5),
    legend.position = "right"
  )

ggsave(file.path(OUTDIR, "HEATMAP_reduced.pdf"), p, width = 10, height = max(4, 0.22 * length(taxa_keep) + 2), device = cairo_pdf)
ggsave(file.path(OUTDIR, "HEATMAP_reduced.png"), p, width = 10, height = max(4, 0.22 * length(taxa_keep) + 2), dpi = 300)

## ---------- 6) 同步导出“精简矩阵”（方便你后续自己做图/做统计） ----------
mat <- plot_df %>%
  select(taxon, metabolite, rho) %>%
  pivot_wider(names_from = metabolite, values_from = rho)

write_csv(mat, file.path(OUTDIR, "HEATMAP_reduced_matrix_rho.csv"))

message("[DONE] Reduced heatmap + tables saved in: ", OUTDIR)
message("Method = ", REDUCE_METHOD, "; taxa kept = ", length(taxa_keep), "; metabolites = ", length(met_set))
