## =========================================================
## STEP_TDC_operon_from_gff3_swissprot.R
## Robust plotting WITHOUT gggenes (avoid unit errors forever)
## - Read GFF3 + SwissProt tab-delimited .xls
## - Find best tdc/tyrDC hit, extract ±WINDOW_BP region
## - Plot gene arrows using ggplot2 polygons (no gggenes)
## - Output PDF (required) + optional PNG + tables
## =========================================================

## ---------------------------
## 0) >>> 可调参数（你主要改这里）<<<
## ---------------------------

## --- 输入文件（工作目录下）---
EFA_GFF <- "Efa13_chr.gff3"
EFA_SW  <- "Efa13_chr_swissprot.xls"

EFI_GFF <- "Efi9_chr.gff3"
EFI_SW  <- "Efi9_chr_swissprot.xls"

## --- 输出文件夹（建议每次改新名字，避免看旧PDF缓存）---
OUTDIR <- "OUT_TDC_operon_polygonArrows_BIG_v1"

## --- operon 窗口大小：以命中中心 ±WINDOW_BP ---
WINDOW_BP <- 10000

## --- 定位 tdc/tyrDC 的关键词（按你 SwissProt 描述可扩充）---
KEY_PAT <- "(\\btyrDC\\b|\\btdcA\\b|\\btdc\\b|tyrosine decarboxylase|aromatic amino acid decarboxylase|\\btyrP\\b|tyramine)"

## --- 标签显示策略 ---
LABEL_MODE <- "smart"
# "smart"   : 标关键基因 + 有 GeneName 的基因（gene_name != "-"）
# "key_only": 只标关键基因（论文更干净，放大字体后推荐）
# "all"     : 窗口内全部标（可能很挤）

KEY_GENES <- c("tyrDC","tdcA","tdc","tyrP","tyrS1","nhaC","yutG","ldh2","pbuG","ydgF")

## --- 全局放大倍数：想“整体大一倍”就 2 ---
SCALE_FACTOR <- 1

## ====== 这三个参数决定“箭头+字体怎么变大” ======
## 1) 箭头高度/粗细（最关键：你要箭头更高就调它）
ARROW_HEIGHT_Y <- min(0.92, 0.42 * SCALE_FACTOR)   # y方向高度（0~1之间），建议 0.6~0.92

## 2) 箭头头部长度（bp）：越大箭头尖越长
ARROW_HEAD_LEN_BP <- round(220 * SCALE_FACTOR)     # 建议 200~600

## 3) 字体大小（ggplot::geom_text 的 size）
LABEL_TEXT_SIZE <- 3 * SCALE_FACTOR              # 基因标签（建议 6~12）
TITLE_TEXT_SIZE <- 10  * SCALE_FACTOR
AXIS_TITLE_SIZE <- 9   * SCALE_FACTOR
AXIS_TEXT_SIZE  <- 8   * SCALE_FACTOR
Y_TEXT_SIZE     <- 9   * SCALE_FACTOR
LABEL_FONTFACE  <- "bold"

## --- 画布大小（英寸）---
PLOT_WIDTH_SINGLE  <- 11
PLOT_HEIGHT_SINGLE <- 2.8 * SCALE_FACTOR/1.15
PLOT_WIDTH_BOTH    <- 11
PLOT_HEIGHT_BOTH   <- 4.2 * SCALE_FACTOR/1.15

## --- 输出格式 ---
MAKE_PNG <- TRUE
DPI_PNG  <- 300

## ---------------------------
## 1) Packages
## ---------------------------
pkg_needed <- c("dplyr","stringr","ggplot2","readr","tibble")
pkg_missing <- pkg_needed[!sapply(pkg_needed, requireNamespace, quietly = TRUE)]
if (length(pkg_missing) > 0) {
  install.packages(pkg_missing, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(readr)
  library(tibble)
})

## ---------------------------
## 2) Checks + output dir
## ---------------------------
stopifnot(file.exists(EFA_GFF), file.exists(EFA_SW), file.exists(EFI_GFF), file.exists(EFI_SW))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

message("=== Size parameters ===")
message("ARROW_HEIGHT_Y = ", ARROW_HEIGHT_Y,
        " | HEAD_LEN_BP = ", ARROW_HEAD_LEN_BP,
        " | LABEL_TEXT_SIZE = ", LABEL_TEXT_SIZE)
message("Output folder: ", normalizePath(OUTDIR))

## ---------------------------
## 3) I/O helpers
## ---------------------------
find_col <- function(df, patterns) {
  nms <- names(df)
  for (pat in patterns) {
    hit <- nms[str_detect(tolower(nms), tolower(pat))]
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

safe_as_numeric <- function(x) suppressWarnings(as.numeric(x))

read_swissprot_table <- function(sw_file, strain) {
  ## SwissProt 结果通常是 TAB 分隔文本（虽然扩展名 .xls）
  sw <- read.delim(sw_file, sep = "\t", header = TRUE, comment.char = "#",
                   stringsAsFactors = FALSE, check.names = FALSE)

  col_orf  <- find_col(sw, c("^orf name$", "orf", "query"))
  col_gene <- find_col(sw, c("genename", "gene name"))
  col_desc <- find_col(sw, c("description", "protein", "annotation", "function"))
  col_eval <- find_col(sw, c("e-value", "evalue", "expect"))
  col_bit  <- find_col(sw, c("bitscore", "bit score", "^score$"))

  if (is.na(col_orf)) stop(paste0("[", strain, "] SwissProt table cannot find ORF name column."))

  sw2 <- sw %>%
    mutate(
      strain = strain,
      ORF = .data[[col_orf]],
      orf_num = safe_as_numeric(str_extract(ORF, "\\d+$")),
      gene_name = if (!is.na(col_gene)) .data[[col_gene]] else NA_character_,
      description = if (!is.na(col_desc)) .data[[col_desc]] else NA_character_,
      evalue = if (!is.na(col_eval)) safe_as_numeric(.data[[col_eval]]) else NA_real_,
      bitscore = if (!is.na(col_bit)) safe_as_numeric(.data[[col_bit]]) else NA_real_
    ) %>%
    filter(!is.na(orf_num))

  ## 每个 ORF 取最佳命中：bitscore 高优先，其次 evalue 低
  if (!all(is.na(sw2$bitscore))) {
    sw_best <- sw2 %>%
      arrange(strain, orf_num, desc(bitscore), evalue) %>%
      group_by(strain, orf_num) %>% slice(1) %>% ungroup()
  } else if (!all(is.na(sw2$evalue))) {
    sw_best <- sw2 %>%
      arrange(strain, orf_num, evalue) %>%
      group_by(strain, orf_num) %>% slice(1) %>% ungroup()
  } else {
    sw_best <- sw2 %>%
      group_by(strain, orf_num) %>% slice(1) %>% ungroup()
  }

  sw_best %>% select(strain, orf_num, ORF, gene_name, description, evalue, bitscore)
}

read_gff_cds <- function(gff_file, strain) {
  gff <- readr::read_tsv(gff_file, comment = "#", col_names = FALSE,
                         show_col_types = FALSE, progress = FALSE)
  if (ncol(gff) < 9) stop(paste0("[", strain, "] GFF3 appears to have <9 columns."))

  gff <- gff[,1:9]
  names(gff) <- c("seqid","source","type","start","end","score","strand","phase","attributes")

  gff %>%
    filter(type == "CDS") %>%
    mutate(
      strain = strain,
      contig = as.character(seqid),
      start  = as.integer(start),
      end    = as.integer(end),
      strand = as.character(strand),
      direction = ifelse(strand == "-", -1, 1),
      ID = str_extract(attributes, "(?<=ID=)[^;]+"),
      Parent = str_extract(attributes, "(?<=Parent=)[^;]+"),
      gene_num = suppressWarnings(as.integer(str_extract(Parent, "(?<=gene)\\d+"))),
      gene_num = ifelse(is.na(gene_num), suppressWarnings(as.integer(str_extract(ID, "\\d+$"))), gene_num)
    ) %>%
    filter(!is.na(gene_num))
}

pick_best_hit <- function(hits) {
  if (!all(is.na(hits$bitscore))) {
    hits %>% arrange(desc(bitscore), evalue) %>% slice(1)
  } else if (!all(is.na(hits$evalue))) {
    hits %>% arrange(evalue) %>% slice(1)
  } else {
    hits %>% slice(1)
  }
}

make_operon_dataset <- function(gff_file, sw_file, strain, window_bp, key_pat, label_mode, key_genes) {
  sw_best <- read_swissprot_table(sw_file, strain)

  cds <- read_gff_cds(gff_file, strain) %>%
    left_join(sw_best, by = c("strain"="strain", "gene_num"="orf_num")) %>%
    mutate(
      ann_text = paste(coalesce(gene_name,""), coalesce(description,""), sep=" | "),
      gene_label = case_when(
        !is.na(gene_name) & gene_name != "" & gene_name != "-" ~ gene_name,
        TRUE ~ paste0("gene", gene_num)
      )
    )

  hits <- cds %>% filter(str_detect(ann_text, regex(key_pat, ignore_case = TRUE)))
  if (nrow(hits) == 0) stop(paste0("[", strain, "] Cannot find tdc/tyrDC hit using KEY_PAT."))

  hits_best <- pick_best_hit(hits)

  center <- floor((hits_best$start + hits_best$end)/2)
  win_start <- max(1, center - window_bp)
  win_end   <- center + window_bp

  win_df <- cds %>%
    filter(contig == hits_best$contig, start <= win_end, end >= win_start) %>%
    mutate(
      xmin = pmax(start, win_start),
      xmax = pmin(end, win_end),
      xmin_rel = xmin - win_start + 1,
      xmax_rel = xmax - win_start + 1,
      is_key_by_pat  = str_detect(ann_text, regex(key_pat, ignore_case = TRUE)),
      is_key_by_list = gene_label %in% key_genes
    ) %>%
    mutate(
      label_show = case_when(
        label_mode == "all" ~ gene_label,
        label_mode == "key_only" & (is_key_by_pat | is_key_by_list) ~ gene_label,
        label_mode == "smart" & (is_key_by_pat | is_key_by_list | (!is.na(gene_name) & gene_name != "" & gene_name != "-")) ~ gene_label,
        TRUE ~ ""
      )
    ) %>%
    arrange(xmin)

  summary_row <- tibble(
    strain = strain,
    contig = hits_best$contig,
    tdc_gene_num = hits_best$gene_num,
    tdc_gene_name = hits_best$gene_label,
    tdc_description = hits_best$description,
    tdc_start = hits_best$start,
    tdc_end   = hits_best$end,
    window_start = win_start,
    window_end   = win_end,
    window_bp    = window_bp
  )

  list(all_cds = cds, tdc_hit = hits_best, win_df = win_df, summary = summary_row)
}

write_tables <- function(obj, strain, outdir) {
  write.csv(obj$tdc_hit, file.path(outdir, paste0(strain, "_tdc_hit_best.csv")),
            row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(obj$win_df, file.path(outdir, paste0(strain, "_window_CDS_table.csv")),
            row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(obj$all_cds, file.path(outdir, paste0(strain, "_all_CDS_with_swissprot.csv")),
            row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(obj$summary, file.path(outdir, paste0(strain, "_window_summary.csv")),
            row.names = FALSE, fileEncoding = "UTF-8")
}

## ---------------------------
## 4) Polygon arrow builder (no gggenes)
## ---------------------------
make_arrow_poly_one <- function(xmin, xmax, direction, y, arrow_height_y, head_len_bp, gene_id, gene_label, label_show) {
  len <- xmax - xmin
  head <- min(head_len_bp, max(30, len * 0.35))  # 头部不至于太短/太长
  h <- arrow_height_y/2

  if (direction >= 0) {
    xb <- max(xmin, xmax - head)
    # (xmin,y-h) -> (xb,y-h) -> (xmax,y) -> (xb,y+h) -> (xmin,y+h)
    df <- tibble(
      gene_id = gene_id,
      gene_label = gene_label,
      label_show = label_show,
      x = c(xmin, xb,  xmax, xb,  xmin),
      y = c(y-h,  y-h, y,   y+h, y+h)
    )
  } else {
    xb <- min(xmax, xmin + head)
    # (xmax,y-h) -> (xb,y-h) -> (xmin,y) -> (xb,y+h) -> (xmax,y+h)
    df <- tibble(
      gene_id = gene_id,
      gene_label = gene_label,
      label_show = label_show,
      x = c(xmax, xb,  xmin, xb,  xmax),
      y = c(y-h,  y-h, y,   y+h, y+h)
    )
  }
  df
}

make_arrow_polys <- function(df, x1_col, x2_col, y_value, arrow_height_y, head_len_bp) {
  df2 <- df %>% mutate(gene_id = row_number())
  polys <- lapply(seq_len(nrow(df2)), function(i) {
    make_arrow_poly_one(
      xmin = df2[[x1_col]][i],
      xmax = df2[[x2_col]][i],
      direction = df2$direction[i],
      y = y_value,
      arrow_height_y = arrow_height_y,
      head_len_bp = head_len_bp,
      gene_id = df2$gene_id[i],
      gene_label = df2$gene_label[i],
      label_show = df2$label_show[i]
    )
  })
  bind_rows(polys)
}

## ---------------------------
## 5) Plot functions
## ---------------------------
base_theme_big <- function() {
  theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.title   = element_text(size = TITLE_TEXT_SIZE, face = "bold"),
      axis.title.x = element_text(size = AXIS_TITLE_SIZE),
      axis.text.x  = element_text(size = AXIS_TEXT_SIZE),
      axis.text.y  = element_text(size = Y_TEXT_SIZE)
    )
}

plot_single <- function(win_df, strain, outprefix, outdir, use_relative = FALSE) {
  y0 <- 1

  if (!use_relative) {
    poly <- make_arrow_polys(win_df, "xmin", "xmax", y0, ARROW_HEIGHT_Y, ARROW_HEAD_LEN_BP)
    label_df <- win_df %>%
      mutate(y = y0, x = (xmin + xmax)/2) %>%
      filter(label_show != "")
    p <- ggplot(poly, aes(x=x, y=y, group=gene_id, fill=gene_label)) +
      geom_polygon(linewidth = 0.2) +
      geom_text(data=label_df, aes(x=x, y=y, label=label_show),
                inherit.aes = FALSE, size = LABEL_TEXT_SIZE, fontface = LABEL_FONTFACE) +
      scale_y_continuous(breaks = y0, labels = strain, limits = c(0.5, 1.5)) +
      labs(x="Genomic position (bp)", y=NULL,
           title=paste0(strain, " putative tdc/tyrDC region (±", WINDOW_BP/1000, " kb)")) +
      base_theme_big()
  } else {
    poly <- make_arrow_polys(win_df, "xmin_rel", "xmax_rel", y0, ARROW_HEIGHT_Y, ARROW_HEAD_LEN_BP)
    label_df <- win_df %>%
      mutate(y = y0, x = (xmin_rel + xmax_rel)/2) %>%
      filter(label_show != "")
    p <- ggplot(poly, aes(x=x, y=y, group=gene_id, fill=gene_label)) +
      geom_polygon(linewidth = 0.2) +
      geom_text(data=label_df, aes(x=x, y=y, label=label_show),
                inherit.aes = FALSE, size = LABEL_TEXT_SIZE, fontface = LABEL_FONTFACE) +
      scale_y_continuous(breaks = y0, labels = strain, limits = c(0.5, 1.5)) +
      labs(x=paste0("Relative position within window (bp) [", 2*WINDOW_BP, " bp window]"),
           y=NULL,
           title=paste0(strain, " putative tdc/tyrDC region (relative coordinates)")) +
      base_theme_big()
  }

  ggsave(file.path(outdir, paste0(outprefix, ".pdf")), p,
         width = PLOT_WIDTH_SINGLE, height = PLOT_HEIGHT_SINGLE, device="pdf", useDingbats=FALSE)
  if (isTRUE(MAKE_PNG)) {
    ggsave(file.path(outdir, paste0(outprefix, ".png")), p,
           width = PLOT_WIDTH_SINGLE, height = PLOT_HEIGHT_SINGLE, dpi = DPI_PNG)
  }
}

plot_both <- function(win_df1, win_df2, outprefix, outdir) {
  ymap <- c(Efa13=1, Efi9=2)

  poly1 <- make_arrow_polys(win_df1, "xmin_rel", "xmax_rel", ymap["Efa13"], ARROW_HEIGHT_Y, ARROW_HEAD_LEN_BP) %>%
    mutate(strain="Efa13")
  poly2 <- make_arrow_polys(win_df2, "xmin_rel", "xmax_rel", ymap["Efi9"],  ARROW_HEIGHT_Y, ARROW_HEAD_LEN_BP) %>%
    mutate(strain="Efi9")
  poly <- bind_rows(poly1, poly2)

  label1 <- win_df1 %>% mutate(strain="Efa13", y=ymap["Efa13"], x=(xmin_rel+xmax_rel)/2) %>% filter(label_show!="")
  label2 <- win_df2 %>% mutate(strain="Efi9",  y=ymap["Efi9"],  x=(xmin_rel+xmax_rel)/2) %>% filter(label_show!="")
  label_df <- bind_rows(label1, label2)

  p <- ggplot(poly, aes(x=x, y=y, group=interaction(strain,gene_id), fill=gene_label)) +
    geom_polygon(linewidth = 0.2) +
    geom_text(data=label_df, aes(x=x, y=y, label=label_show),
              inherit.aes = FALSE, size = LABEL_TEXT_SIZE, fontface = LABEL_FONTFACE) +
    scale_y_continuous(breaks = c(1,2), labels = c("Efa13","Efi9"), limits = c(0.5, 2.5)) +
    labs(x=paste0("Relative position within window (bp) [", 2*WINDOW_BP, " bp window]"),
         y=NULL,
         title="Putative tdc/tyrDC regions (Efa13 vs Efi9; aligned by relative coordinates)") +
    base_theme_big()

  ggsave(file.path(outdir, paste0(outprefix, ".pdf")), p,
         width = PLOT_WIDTH_BOTH, height = PLOT_HEIGHT_BOTH, device="pdf", useDingbats=FALSE)
  if (isTRUE(MAKE_PNG)) {
    ggsave(file.path(outdir, paste0(outprefix, ".png")), p,
           width = PLOT_WIDTH_BOTH, height = PLOT_HEIGHT_BOTH, dpi = DPI_PNG)
  }
}

## ---------------------------
## 6) Run
## ---------------------------
efa <- make_operon_dataset(EFA_GFF, EFA_SW, "Efa13", WINDOW_BP, KEY_PAT, LABEL_MODE, KEY_GENES)
efi <- make_operon_dataset(EFI_GFF, EFI_SW, "Efi9",  WINDOW_BP, KEY_PAT, LABEL_MODE, KEY_GENES)

write_tables(efa, "Efa13", OUTDIR)
write_tables(efi, "Efi9",  OUTDIR)

write.csv(bind_rows(efa$summary, efi$summary), file.path(OUTDIR, "BOTH_window_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(bind_rows(efa$win_df, efi$win_df), file.path(OUTDIR, "BOTH_window_CDS_table.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

plot_single(efa$win_df, "Efa13", "Efa13_tdc_operon_window_BIG", OUTDIR, use_relative = FALSE)
plot_single(efi$win_df, "Efi9",  "Efi9_tdc_operon_window_BIG",  OUTDIR, use_relative = FALSE)

plot_both(
  efa$win_df %>% mutate(strain="Efa13"),
  efi$win_df %>% mutate(strain="Efi9"),
  "BOTH_tdc_operon_window_aligned_relative_BIG",
  OUTDIR
)

sink(file.path(OUTDIR, "sessionInfo.txt"))
print(sessionInfo())
sink()

message("Done. Outputs saved in: ", normalizePath(OUTDIR))
