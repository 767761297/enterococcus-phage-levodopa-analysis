## STEP_forest_signature_6Levels_GR_vs_MR_q0.1.R
## Forest plots for fingerprint taxa at 6 taxonomic levels (GR vs MR)
## Effect size = mean(CLR_log2 in GR) - mean(CLR_log2 in MR)
## CI = bootstrap 95% CI within groups
## q = BH-FDR from Wilcoxon tests across ALL taxa at that level

options(error = NULL)  # 避免递归 error handler 引发“严重 bug”

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

# ---------------- USER SETTINGS ----------------
ABUND_FILE <- "allsample.all.reabundance.xls"     # 你的输入文件（有 Phylum~Species 六列 + 样本列）
SIG_DIR    <- "4. 16S_Signature_q0.1_6Levels"        # （可选）指纹集合所在文件夹
USE_SIGNATURE_FILES <- TRUE                        # TRUE: 优先用 SIGNATURE 文件；FALSE: 直接按 q<0.1 定义指纹

Q_CUT   <- 0.10
BOOT_B  <- 800     # bootstrap次数；论文图建议 800~2000；你电脑慢可先 300 试跑
SEED    <- 20260106

OUTDIR <- "16S_Forest_Signature_GR_vs_MR_q0.1_6Levels"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

LEVELS <- c("Phylum","Class","Order","Family","Genus","Species")

# 颜色
COL_UP   <- "#D7191C"   # GR > MR
COL_DOWN <- "#2C7BB6"   # GR < MR
# ----------------------------------------------


# ---------- helpers ----------
log_msg <- function(...) cat(paste0(..., "\n"))

read_tsv_loose <- function(path) {
  # 兼容 .xls(实际为tab分隔)；避免编码问题
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
             stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
}

infer_group <- function(sample_ids) {
  # 默认按首字母：A=GR, B=MR, C=MR_Phage, D=MR_AFMT
  prefix <- substring(sample_ids, 1, 1)
  grp <- dplyr::case_when(
    prefix == "A" ~ "GR",
    prefix == "B" ~ "MR",
    prefix == "C" ~ "MR+Phage cocktail",
    prefix == "D" ~ "MR+AFMT",
    TRUE ~ NA_character_
  )
  if (any(is.na(grp))) {
    bad <- unique(sample_ids[is.na(grp)])
    stop("Cannot infer group for samples: ", paste(bad, collapse = ", "),
         "\nPlease rename samples with A/B/C/D prefix or provide your own mapping in infer_group().")
  }
  grp
}

clr_log2 <- function(mat, pseudo = NULL) {
  # mat: taxa x samples, non-negative
  x <- as.matrix(mat)
  x[x < 0] <- 0
  nz <- x[x > 0]
  if (is.null(pseudo)) {
    if (length(nz) == 0) pseudo <- 1e-6 else pseudo <- max(min(nz) / 2, 1e-12)
  }
  lx <- log(x + pseudo)
  # CLR per sample: subtract mean log across taxa
  lx_center <- sweep(lx, 2, colMeans(lx), FUN = "-")
  lx_center / log(2)
}

wilcox_p <- function(x, y) {
  # 稳健：遇到全相等/常量返回 1
  out <- tryCatch({
    stats::wilcox.test(x, y, alternative = "two.sided", exact = FALSE)$p.value
  }, error = function(e) 1)
  if (is.na(out)) out <- 1
  out
}

bootstrap_ci <- function(x_gr, x_mr, B = 800, seed = 1) {
  set.seed(seed)
  ng <- length(x_gr); nm <- length(x_mr)
  if (ng < 2 || nm < 2) return(c(NA_real_, NA_real_))
  boot <- replicate(B, {
    mean(sample(x_gr, ng, replace = TRUE)) - mean(sample(x_mr, nm, replace = TRUE))
  })
  stats::quantile(boot, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
}

read_signature_taxa <- function(level) {
  f <- file.path(SIG_DIR, paste0("SIGNATURE_q0.1_", tolower(level), "_GR_vs_MR.csv"))
  if (!file.exists(f)) return(character(0))
  sig <- tryCatch(read.csv(f, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.null(sig) || !("taxon" %in% names(sig))) return(character(0))
  unique(as.character(sig$taxon))
}

# ---------- main ----------
stopifnot(file.exists(ABUND_FILE))
df0 <- read_tsv_loose(ABUND_FILE)

need_cols <- LEVELS
if (!all(need_cols %in% names(df0))) {
  stop("Input file must contain columns: ", paste(need_cols, collapse = ", "),
       "\nDetected columns: ", paste(names(df0), collapse = ", "))
}

sample_cols <- setdiff(names(df0), LEVELS)
if (length(sample_cols) < 4) stop("No sample columns detected in input file.")

meta <- data.frame(
  sample = sample_cols,
  group  = infer_group(sample_cols),
  stringsAsFactors = FALSE
)

# 只做 GR vs MR
meta2 <- meta %>% filter(group %in% c("GR","MR"))
GR_samps <- meta2$sample[meta2$group == "GR"]
MR_samps <- meta2$sample[meta2$group == "MR"]
if (length(GR_samps) < 2 || length(MR_samps) < 2) stop("Need >=2 samples in both GR and MR.")

log_msg("[INFO] Samples: GR=", length(GR_samps), " | MR=", length(MR_samps))
log_msg("[INFO] Output dir: ", OUTDIR)

all_tables <- list()

for (lv in LEVELS) {
  log_msg("\n[LEVEL] ", lv, " ...")

  # 聚合到该层级
  df_lv <- df0 %>%
    transmute(
      taxon = as.character(.data[[lv]]),
      across(all_of(sample_cols), ~ suppressWarnings(as.numeric(.x)))
    ) %>%
    filter(!is.na(taxon) & taxon != "") %>%
    group_by(taxon) %>%
    summarise(across(all_of(sample_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

  mat <- as.matrix(df_lv[, sample_cols, drop = FALSE])
  rownames(mat) <- df_lv$taxon

  # 仅保留 GR/MR 样本列
  mat2 <- mat[, c(GR_samps, MR_samps), drop = FALSE]

  # 去掉在 GR+MR 全为0的taxa
  keep <- rowSums(mat2) > 0
  mat2 <- mat2[keep, , drop = FALSE]
  if (nrow(mat2) < 3) {
    log_msg("[SKIP] too few taxa after filtering at ", lv)
    next
  }

  # CLR (log2)
  clr <- clr_log2(mat2)

  # 统计：Wilcoxon + BH
  taxa <- rownames(clr)
  pvals <- vapply(taxa, function(t) {
    xg <- clr[t, GR_samps]
    xm <- clr[t, MR_samps]
    wilcox_p(xg, xm)
  }, numeric(1))
  qvals <- p.adjust(pvals, method = "BH")

  # 效应量（Δmean CLR_log2）
  eff <- vapply(taxa, function(t) mean(clr[t, GR_samps]) - mean(clr[t, MR_samps]), numeric(1))

  stats_all <- data.frame(
    level = lv,
    taxon = taxa,
    effect = eff,
    p_value = pvals,
    q_value = qvals,
    direction = ifelse(eff >= 0, "GR higher", "MR higher"),
    stringsAsFactors = FALSE
  )

  # 指纹集合：优先读 SIGNATURE 文件；否则 q<0.1
  sig_taxa <- character(0)
  if (USE_SIGNATURE_FILES) sig_taxa <- read_signature_taxa(lv)

  if (length(sig_taxa) > 0) {
    sig <- stats_all %>% filter(taxon %in% sig_taxa)
    log_msg("[INFO] signature from file: n=", nrow(sig))
  } else {
    sig <- stats_all %>% filter(q_value < Q_CUT)
    log_msg("[INFO] signature from q<", Q_CUT, ": n=", nrow(sig))
  }

  if (nrow(sig) == 0) {
    log_msg("[SKIP] no signature taxa at ", lv)
    next
  }

  # bootstrap CI（只对指纹做，避免太慢）
  set.seed(SEED)
  cis <- t(vapply(sig$taxon, function(t) {
    xg <- clr[t, GR_samps]
    xm <- clr[t, MR_samps]
    bootstrap_ci(xg, xm, B = BOOT_B, seed = SEED + match(t, sig$taxon))
  }, numeric(2)))
  sig$ci_low  <- cis[,1]
  sig$ci_high <- cis[,2]

  # 排序：按效应量
  sig <- sig %>%
    arrange(effect) %>%
    mutate(taxon = factor(taxon, levels = taxon))

  # 输出表
  out_csv <- file.path(OUTDIR, paste0("Forest_signature_q0.1_", tolower(lv), "_GR_vs_MR.csv"))
  write.csv(sig, out_csv, row.names = FALSE)

  all_tables[[lv]] <- sig

  # 画森林图（动态高度，避免挤）
  n_tax <- nrow(sig)
  h <- max(4.5, 0.22 * n_tax + 1.2)

  p <- ggplot(sig, aes(x = effect, y = taxon)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey40") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high, color = direction), height = 0.2, linewidth = 0.7) +
    geom_point(aes(color = direction), size = 2.4, shape = 16) +
    scale_color_manual(values = c("GR higher" = COL_UP, "MR higher" = COL_DOWN)) +
    theme_classic(base_size = 12) +
    theme(
      legend.title = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y  = element_text(size = 9)
    ) +
    labs(x = "Effect size")

  out_pdf <- file.path(OUTDIR, paste0("Forest_signature_q0.1_", tolower(lv), "_GR_vs_MR.pdf"))
  out_png <- file.path(OUTDIR, paste0("Forest_signature_q0.1_", tolower(lv), "_GR_vs_MR.png"))

  ggsave(out_pdf, p, width = 8.2, height = h)
  ggsave(out_png, p, width = 8.2, height = h, dpi = 300)

  log_msg("[OK] ", lv, " -> ", basename(out_pdf), " | table: ", basename(out_csv))
}

# 可选：汇总输出（如果安装了 writexl）
if (length(all_tables) > 0 && requireNamespace("writexl", quietly = TRUE)) {
  out_xlsx <- file.path(OUTDIR, "Forest_signature_q0.1_6Levels_GR_vs_MR.xlsx")
  writexl::write_xlsx(all_tables, out_xlsx)
  log_msg("\n[OK] Combined xlsx: ", out_xlsx)
} else {
  log_msg("\n[NOTE] writexl not installed -> skip combined xlsx (CSV per level already written).")
}

log_msg("\n[DONE] ", OUTDIR)
