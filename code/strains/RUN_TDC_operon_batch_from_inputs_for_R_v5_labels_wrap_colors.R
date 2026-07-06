## =========================================================
## RUN_TDC_operon_batch_from_inputs_for_R_v3_beautified.R
## - 更稳健地读取 assembly_data_report.txt，优先显示物种/菌株名而不是 accession
## - 多株总图改为 facet 纵向排布，避免 y 轴拥挤
## - 只标注 tdc/tyrDC 及其关键上下游几个基因，避免满图文字重叠
## - 简化配色：关键基因高亮，其余邻近基因浅灰，其它基因更淡
## - 保留 polygon arrow 风格
## =========================================================

## ---------------------------
## 0) >>> 你主要改这里 <<<
## ---------------------------
INPUTS_DIR <- "C:/Users/76776/Desktop/human_full_lenth_16s_phage project/tyrDC_operon_genomes/inputs_for_R"
REPORT_FILE <- "C:/Users/76776/Desktop/human_full_lenth_16s_phage project/tyrDC_operon_genomes/assembly_data_report.txt"

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTDIR <- file.path(dirname(INPUTS_DIR), paste0("results_tyrDC_operon_batch_v3_beautified_", timestamp))

WINDOW_BP <- 10000  # ±10 kb

## 目标关键词：用于寻找 focal tdc/tyrDC 基因
KEY_PAT <- "(\\btyrDC\\b|\\btdcA\\b|\\btdc\\b|tyrosine decarboxylase|aromatic amino acid decarboxylase|\\btyrP\\b|tyramine)"

## operon 里常见、值得额外高亮/标注的基因
KEY_GENES <- c("tyrDC","tdcA","tdc","tyrP","tyrS","tyrS1","nhaC","pbuG","ydgF","adhB","cla-hy","yutG","ldh2","cls")

## 只标注 focal 基因上下游的几个关键基因
LABEL_NEIGHBORS_EACH_SIDE <- 3   # focal 左右各标几个基因，可改成 2/4/5
LABEL_MAX_NCHAR <- 40            # 标签总长度上限（超过则截断）
LABEL_WRAP_WIDTH <- 16           # 长标签换行宽度
LABEL_STYLE <- "organism_strain" # 可选: "organism_strain" 或 "strain_only"

## 图形参数
SCALE_FACTOR <- 1
ARROW_HEIGHT_Y <- 0.54
ARROW_HEAD_LEN_BP <- round(220 * SCALE_FACTOR)
LABEL_TEXT_SIZE <- 3.0 * SCALE_FACTOR
TITLE_TEXT_SIZE <- 12
AXIS_TITLE_SIZE <- 10
AXIS_TEXT_SIZE  <- 8.5
STRIP_TEXT_SIZE <- 8.8
LABEL_FONTFACE  <- "bold"

PLOT_WIDTH_SINGLE  <- 12
PLOT_HEIGHT_SINGLE <- 3.8
PLOT_WIDTH_MULTI   <- 13
MAKE_PNG <- TRUE
DPI_PNG  <- 320

## ---------------------------
## 1) Packages
## ---------------------------
pkg_needed <- c("dplyr","stringr","ggplot2","readr","tibble","ggrepel","jsonlite")
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
  library(ggrepel)
  library(jsonlite)
})

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "per_accession"), showWarnings = FALSE, recursive = TRUE)

## ---------------------------
## 2) Helpers
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

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

extract_gn <- function(x) {
  m <- str_match(coalesce(x, ""), "GN=([^ ]+)")
  m[,2]
}

get_attr <- function(attr, key) {
  m <- str_match(attr, paste0("(?:^|;)", key, "=([^;]+)"))
  m[,2]
}

truncate_label <- function(x, max_n = LABEL_MAX_NCHAR) {
  x <- coalesce(as.character(x), "")
  x <- str_squish(str_replace_all(x, "\\s+", " "))
  ifelse(nchar(x) > max_n, paste0(substr(x, 1, max_n - 3), "..."), x)
}

wrap_plot_label <- function(x, width = LABEL_WRAP_WIDTH, max_n = LABEL_MAX_NCHAR) {
  x <- truncate_label(x, max_n = max_n)
  x <- str_replace_all(x, "([/_-])", "\\1 ")
  x <- str_wrap(x, width = width)
  x <- str_replace_all(x, "([/_-])\\s+", "\\1")
  x
}

pick_display_label <- function(gff_gene, gene_name, gff_product, gene_num) {
  out <- dplyr::case_when(
    !is.na(gff_gene) & gff_gene != "" ~ gff_gene,
    !is.na(gene_name) & gene_name != "" ~ gene_name,
    !is.na(gff_product) & gff_product != "" ~ gff_product,
    TRUE ~ paste0("gene", gene_num)
  )
  str_squish(str_replace_all(out, "\\s+", " "))
}

make_pretty_label <- function(accession, organism = NA_character_, strain = NA_character_, style = LABEL_STYLE) {
  accession <- as.character(accession)
  organism  <- coalesce(as.character(organism), "")
  strain    <- coalesce(as.character(strain), "")

  organism <- str_squish(organism)
  strain   <- str_squish(str_replace(strain, regex("^strain\\s*=\\s*", ignore_case = TRUE), ""))

  has_org <- organism != ""
  has_strain <- strain != ""

  out <- if (identical(style, "strain_only")) {
    ifelse(has_strain, strain, ifelse(has_org, organism, accession))
  } else {
    ifelse(has_org, organism, ifelse(has_strain, strain, accession))
  }

  out[out == ""] <- accession[out == ""]
  out
}

read_swissprot_table <- function(sw_file, strain_id) {
  sw <- read.delim(sw_file, sep = "\t", header = TRUE, comment.char = "#",
                   stringsAsFactors = FALSE, check.names = FALSE)

  col_orf  <- find_col(sw, c("^orf name$", "orf", "query", "^ORF$"))
  col_desc <- find_col(sw, c("description", "stitle", "title", "protein", "annotation", "function"))
  col_eval <- find_col(sw, c("e-value", "evalue", "expect"))
  col_bit  <- find_col(sw, c("bitscore", "bit score", "^score$"))

  if (is.na(col_orf)) stop(paste0("[", strain_id, "] SwissProt table cannot find ORF/query column."))

  sw2 <- sw %>%
    mutate(
      strain_id = strain_id,
      ORF = .data[[col_orf]],
      orf_num = safe_as_numeric(str_extract(ORF, "\\d+$")),
      description = if (!is.na(col_desc)) .data[[col_desc]] else NA_character_,
      gene_name = extract_gn(description),
      evalue = if (!is.na(col_eval)) safe_as_numeric(.data[[col_eval]]) else NA_real_,
      bitscore = if (!is.na(col_bit)) safe_as_numeric(.data[[col_bit]]) else NA_real_
    ) %>%
    filter(!is.na(orf_num))

  if (!all(is.na(sw2$bitscore))) {
    sw_best <- sw2 %>%
      arrange(strain_id, orf_num, desc(bitscore), evalue) %>%
      group_by(strain_id, orf_num) %>% slice(1) %>% ungroup()
  } else if (!all(is.na(sw2$evalue))) {
    sw_best <- sw2 %>%
      arrange(strain_id, orf_num, evalue) %>%
      group_by(strain_id, orf_num) %>% slice(1) %>% ungroup()
  } else {
    sw_best <- sw2 %>%
      group_by(strain_id, orf_num) %>% slice(1) %>% ungroup()
  }

  sw_best %>% select(strain_id, orf_num, ORF, gene_name, description, evalue, bitscore)
}

read_gff_cds <- function(gff_file, strain_id) {
  gff <- readr::read_tsv(gff_file, comment = "#", col_names = FALSE,
                         show_col_types = FALSE, progress = FALSE)
  if (ncol(gff) < 9) stop(paste0("[", strain_id, "] GFF3 appears to have <9 columns."))
  gff <- gff[,1:9]
  names(gff) <- c("seqid","source","type","start","end","score","strand","phase","attributes")

  gff %>%
    filter(type == "CDS") %>%
    mutate(
      strain_id = strain_id,
      contig = as.character(seqid),
      start  = as.integer(start),
      end    = as.integer(end),
      strand = as.character(strand),
      direction = ifelse(strand == "-", -1, 1),
      ID = get_attr(attributes, "ID"),
      Parent = get_attr(attributes, "Parent"),
      gff_gene = get_attr(attributes, "gene"),
      gff_product = get_attr(attributes, "product"),
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


read_accession_labels <- function(report_file) {
  if (!file.exists(report_file)) {
    message("[WARN] assembly_data_report.txt not found: ", report_file)
    return(tibble(accession=character(), organism=character(), strain=character(), label=character()))
  }

  lines <- readLines(report_file, warn = FALSE, encoding = "UTF-8")
  lines <- sub("^\ufeff", "", lines)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0) {
    return(tibble(accession=character(), organism=character(), strain=character(), label=character()))
  }

  ## 情况 1：JSONL（你这个 assembly_data_report.txt 就是这种）
  if (grepl("^\\{", trimws(lines[1]))) {
    parse_one_json_line <- function(ln) {
      obj <- tryCatch(jsonlite::fromJSON(ln, simplifyVector = FALSE), error = function(e) NULL)
      if (is.null(obj)) return(NULL)

      acc <- obj$currentAccession %||% obj$accession %||% obj$pairedAccession
      if (is.null(acc) && !is.null(obj$pairedAssembly$accession)) acc <- obj$pairedAssembly$accession

      organism <- obj$organism$organismName %||%
        obj$biosample$description$organism$organismName %||%
        obj$averageNucleotideIdentity$submittedOrganism %||% NA_character_

      strain <- obj$organism$infraspecificNames$strain %||%
        obj$biosample$strain %||% NA_character_

      if ((is.null(strain) || is.na(strain) || strain == "") && !is.null(obj$biosample$attributes)) {
        attrs <- obj$biosample$attributes
        if (is.list(attrs) && length(attrs) > 0) {
          idx <- which(vapply(attrs, function(z) identical(tolower(z$name %||% ""), "strain"), logical(1)))
          if (length(idx) > 0) {
            strain <- attrs[[idx[1]]]$value %||% strain
          }
        }
      }

      tibble(
        accession = as.character(acc %||% NA_character_),
        organism  = as.character(organism %||% NA_character_),
        strain    = as.character(strain %||% NA_character_)
      )
    }

    df <- lapply(lines, parse_one_json_line) %>%
      bind_rows() %>%
      filter(!is.na(accession), accession != "") %>%
      distinct(accession, .keep_all = TRUE)

    return(df %>% mutate(label = make_pretty_label(accession, organism, strain)))
  }

  ## 情况 2：标准 NCBI assembly_data_report.txt（注释行 + 表头行）
  comment_idx <- grep("^#", lines)
  parsed <- NULL

  if (length(comment_idx) > 0) {
    header_line <- lines[max(comment_idx)]
    header_line <- sub("^#\\s*", "", header_line)
    header_fields <- strsplit(header_line, "\t", fixed = TRUE)[[1]]

    data_lines <- lines[-seq_len(max(comment_idx))]
    data_lines <- data_lines[nzchar(trimws(data_lines))]

    if (length(header_fields) > 1 && length(data_lines) > 0) {
      tmp <- tempfile(fileext = ".tsv")
      writeLines(c(paste(header_fields, collapse = "\t"), data_lines), tmp, useBytes = TRUE)
      parsed <- tryCatch(readr::read_tsv(tmp, show_col_types = FALSE), error = function(e) NULL)
    }
  }

  ## 情况 3：普通 TSV
  if (is.null(parsed)) {
    parsed <- tryCatch(readr::read_tsv(report_file, show_col_types = FALSE, comment = ""), error = function(e) NULL)
  }

  if (is.null(parsed) || nrow(parsed) == 0) {
    message("[WARN] Cannot parse assembly_data_report.txt; y labels will fall back to accession.")
    return(tibble(accession=character(), organism=character(), strain=character(), label=character()))
  }

  cn <- names(parsed)
  cn_low <- tolower(cn)

  acc_idx <- match(TRUE, cn_low %in% c("accession", "assembly accession", "assembly_accession", "assembly"), nomatch = 0)
  acc_col <- if (acc_idx > 0) cn[acc_idx] else character(0)
  if (length(acc_col) == 0) {
    acc_hits <- which(str_detect(cn_low, "accession"))
    acc_col <- if (length(acc_hits) > 0) cn[acc_hits[1]] else character(0)
  }

  org_hits1 <- which(cn_low %in% c("organism name", "organism_name", "organism"))
  org_col <- if (length(org_hits1) > 0) cn[org_hits1[1]] else character(0)
  if (length(org_col) == 0) {
    org_hits2 <- which(str_detect(cn_low, "organism"))
    org_col <- if (length(org_hits2) > 0) cn[org_hits2[1]] else character(0)
  }

  strain_hits1 <- which(cn_low %in% c("infraspecific name", "infraspecific_name", "strain", "isolate"))
  strain_col <- if (length(strain_hits1) > 0) cn[strain_hits1[1]] else character(0)
  if (length(strain_col) == 0) {
    strain_hits2 <- which(str_detect(cn_low, "strain|infraspecific|isolate"))
    strain_col <- if (length(strain_hits2) > 0) cn[strain_hits2[1]] else character(0)
  }

  if (length(acc_col) == 0) {
    message("[WARN] accession column not found in assembly_data_report.txt; y labels will fall back to accession.")
    return(tibble(accession=character(), organism=character(), strain=character(), label=character()))
  }

  tibble(
    accession = as.character(parsed[[acc_col]]),
    organism  = if (length(org_col) > 0) as.character(parsed[[org_col]]) else NA_character_,
    strain    = if (length(strain_col) > 0) as.character(parsed[[strain_col]]) else NA_character_
  ) %>%
    filter(!is.na(accession), accession != "") %>%
    distinct(accession, .keep_all = TRUE) %>%
    mutate(label = make_pretty_label(accession, organism, strain))
}

make_operon_dataset <- function(gff_file, sw_file, strain_id, window_bp, key_pat, key_genes,
                                label_neighbors_each_side = 3) {
  sw_best <- read_swissprot_table(sw_file, strain_id)

  cds <- read_gff_cds(gff_file, strain_id) %>%
    left_join(sw_best, by = c("strain_id" = "strain_id", "gene_num" = "orf_num")) %>%
    mutate(
      ann_text = paste(coalesce(gff_gene, ""), coalesce(gff_product, ""),
                       coalesce(gene_name, ""), coalesce(description, ""), sep = " | "),
      gene_label = pick_display_label(gff_gene, gene_name, gff_product, gene_num)
    )

  hits <- cds %>% filter(str_detect(ann_text, regex(key_pat, ignore_case = TRUE)))
  if (nrow(hits) == 0) stop(paste0("[", strain_id, "] Cannot find tdc/tyrDC hit using KEY_PAT."))

  hits_best <- pick_best_hit(hits)
  center <- floor((hits_best$start + hits_best$end) / 2)
  win_start <- max(1, center - window_bp)
  win_end   <- center + window_bp

  win_df <- cds %>%
    filter(contig == hits_best$contig, start <= win_end, end >= win_start) %>%
    mutate(
      xmin = pmax(start, win_start),
      xmax = pmin(end, win_end),
      xmin_rel = xmin - win_start + 1,
      xmax_rel = xmax - win_start + 1,
      gene_center = (xmin + xmax) / 2,
      gene_center_rel = (xmin_rel + xmax_rel) / 2,
      is_hit_by_pat = str_detect(ann_text, regex(key_pat, ignore_case = TRUE)),
      is_key_gene = gene_label %in% key_genes | str_detect(ann_text, regex("\\b(tyrDC|tdcA|tdc|tyrP|tyrS|nhaC)\\b", ignore_case = TRUE))
    ) %>%
    arrange(xmin) %>%
    mutate(row_idx = row_number())

  focal_idx <- which(win_df$gene_num == hits_best$gene_num)
  if (length(focal_idx) == 0) {
    focal_idx <- which.min(abs(win_df$gene_center - center))
  }
  focal_idx <- focal_idx[1]

  neighbor_idx <- seq(max(1, focal_idx - label_neighbors_each_side),
                      min(nrow(win_df), focal_idx + label_neighbors_each_side))
  ## 只保留 focal 及其局部上下游几个基因，避免整窗所有 key gene 都被标出来
  local_key_idx <- which(win_df$is_key_gene & abs(win_df$row_idx - focal_idx) <= (label_neighbors_each_side + 1))
  label_idx <- sort(unique(c(focal_idx, neighbor_idx, local_key_idx)))

  win_df <- win_df %>%
    mutate(
      is_focal = row_idx == focal_idx,
      is_label_target = row_idx %in% label_idx,
      label_show = ifelse(is_label_target, wrap_plot_label(gene_label), ""),
      fill_group = gene_label
    )
  summary_row <- tibble(
    strain_id = strain_id,
    contig = hits_best$contig,
    tdc_gene_num = hits_best$gene_num,
    tdc_gene_name = hits_best$gene_label,
    tdc_description = hits_best$description,
    tdc_start = hits_best$start,
    tdc_end   = hits_best$end,
    focal_center_bp = center,
    window_start = win_start,
    window_end   = win_end,
    window_bp    = window_bp
  )

  list(all_cds = cds, tdc_hit = hits_best, win_df = win_df, summary = summary_row)
}

## Polygon arrow builder
make_arrow_poly_one <- function(xmin, xmax, direction, y, arrow_height_y, head_len_bp,
                                gene_id, gene_label, label_show, fill_group, is_focal) {
  len <- xmax - xmin
  head <- min(head_len_bp, max(30, len * 0.35))
  h <- arrow_height_y / 2

  if (direction >= 0) {
    xb <- max(xmin, xmax - head)
    tibble(
      gene_id = gene_id,
      gene_label = gene_label,
      label_show = label_show,
      fill_group = fill_group,
      is_focal = is_focal,
      x = c(xmin, xb, xmax, xb, xmin),
      y = c(y-h, y-h, y, y+h, y+h)
    )
  } else {
    xb <- min(xmax, xmin + head)
    tibble(
      gene_id = gene_id,
      gene_label = gene_label,
      label_show = label_show,
      fill_group = fill_group,
      is_focal = is_focal,
      x = c(xmax, xb, xmin, xb, xmax),
      y = c(y-h, y-h, y, y+h, y+h)
    )
  }
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
      label_show = df2$label_show[i],
      fill_group = df2$fill_group[i],
      is_focal = df2$is_focal[i]
    )
  })
  bind_rows(polys)
}

base_theme_single <- function() {
  theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = TITLE_TEXT_SIZE, face = "bold"),
      axis.title.x = element_text(size = AXIS_TITLE_SIZE),
      axis.text.x = element_text(size = AXIS_TEXT_SIZE),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(12, 18, 10, 18)
    )
}

base_theme_multi <- function() {
  theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.y = grid::unit(0.25, "lines"),
      strip.placement = "outside",
      strip.background = element_rect(fill = "grey97", color = "grey85"),
      strip.text.y.left = element_text(size = STRIP_TEXT_SIZE, face = "bold", angle = 0),
      plot.title = element_text(size = TITLE_TEXT_SIZE, face = "bold"),
      axis.title.x = element_text(size = AXIS_TITLE_SIZE),
      axis.text.x = element_text(size = AXIS_TEXT_SIZE),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(12, 18, 10, 18)
    )
}

build_label_df <- function(df, x_col) {
  df %>%
    filter(label_show != "") %>%
    mutate(
      x = .data[[x_col]],
      y = 1,
      label_show = as.character(label_show)
    )
}

plot_single <- function(win_df, summary_row, accession, pretty_label, outprefix, outdir) {
  poly <- make_arrow_polys(win_df, "xmin", "xmax", 1, ARROW_HEIGHT_Y, ARROW_HEAD_LEN_BP)
  lab_df <- build_label_df(win_df %>% mutate(x_anchor = gene_center), "x_anchor")

  p <- ggplot(poly, aes(x = x, y = y, group = gene_id, fill = fill_group)) +
    geom_polygon(color = "white", linewidth = 0.25) +
    geom_polygon(data = poly %>% filter(is_focal), fill = NA, color = "black", linewidth = 0.55) +
    geom_vline(xintercept = summary_row$focal_center_bp, linetype = "dashed", color = "grey55", linewidth = 0.35) +
    geom_text_repel(
      data = lab_df,
      aes(x = x, y = y, label = label_show),
      inherit.aes = FALSE,
      direction = "x",
      nudge_y = 0.28,
      box.padding = 0.22,
      point.padding = 0.08,
      segment.color = "grey55",
      segment.size = 0.28,
      min.segment.length = 0,
      max.overlaps = Inf,
      force = 0.8,
      size = LABEL_TEXT_SIZE,
      fontface = LABEL_FONTFACE,
      seed = 1
    ) +
    coord_cartesian(ylim = c(0.58, 1.72), clip = "off") +
    labs(
      x = "Genomic position (bp)",
      y = NULL,
      title = paste0(pretty_label, " (", accession, ")  putative tdc/tyrDC region (±", WINDOW_BP/1000, " kb)")
    ) +
    base_theme_single()

  ggsave(file.path(outdir, paste0(outprefix, ".pdf")), p,
         width = PLOT_WIDTH_SINGLE, height = PLOT_HEIGHT_SINGLE,
         device = "pdf", useDingbats = FALSE)
  if (isTRUE(MAKE_PNG)) {
    ggsave(file.path(outdir, paste0(outprefix, ".png")), p,
           width = PLOT_WIDTH_SINGLE, height = PLOT_HEIGHT_SINGLE, dpi = DPI_PNG)
  }
}

## ---------------------------
## 3) 批量运行
## ---------------------------
acc_dirs <- list.dirs(INPUTS_DIR, recursive = FALSE, full.names = TRUE)
acc_dirs <- acc_dirs[grepl("(GCF|GCA)_\\d+\\.\\d+$", basename(acc_dirs))]
if (length(acc_dirs) == 0) stop("No accession folders found under INPUTS_DIR.")

labels_map <- read_accession_labels(REPORT_FILE)

manifest <- tibble(
  accession = basename(acc_dirs),
  dir = acc_dirs,
  gff = file.path(acc_dirs, paste0(basename(acc_dirs), "_chr.gff3")),
  sw  = file.path(acc_dirs, paste0(basename(acc_dirs), "_chr_swissprot.xls"))
) %>%
  left_join(labels_map %>% select(accession, organism, strain, label), by = "accession") %>%
  mutate(label = ifelse(is.na(label) | label == "", accession, label))

manifest_ok <- manifest %>% filter(file.exists(gff) & file.exists(sw))
if (nrow(manifest_ok) == 0) stop("No valid (gff3 + swissprot) pairs found.")

## 先按物种/菌株标签排序；如无标签则按 accession
manifest_ok <- manifest_ok %>% arrange(label, accession)

write.csv(manifest_ok, file.path(OUTDIR, "manifest_used.csv"), row.names = FALSE, fileEncoding = "UTF-8")

all_summaries <- list()
all_windows <- list()

for (i in seq_len(nrow(manifest_ok))) {
  acc <- manifest_ok$accession[i]
  lab <- manifest_ok$label[i]
  gff <- manifest_ok$gff[i]
  sw  <- manifest_ok$sw[i]

  message("[RUN] ", acc, " | ", str_replace_all(lab, "\\n", " "))

  out_acc <- file.path(OUTDIR, "per_accession", acc)
  dir.create(out_acc, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(out_acc, "tables"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(out_acc, "figures"), showWarnings = FALSE, recursive = TRUE)

  obj <- make_operon_dataset(
    gff_file = gff,
    sw_file = sw,
    strain_id = acc,
    window_bp = WINDOW_BP,
    key_pat = KEY_PAT,
    key_genes = KEY_GENES,
    label_neighbors_each_side = LABEL_NEIGHBORS_EACH_SIDE
  )

  write.csv(obj$tdc_hit, file.path(out_acc, "tables", paste0(acc, "_tdc_hit_best.csv")),
            row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(obj$win_df, file.path(out_acc, "tables", paste0(acc, "_window_CDS_table.csv")),
            row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(obj$summary %>% mutate(accession = acc, label = lab),
            file.path(out_acc, "tables", paste0(acc, "_window_summary.csv")),
            row.names = FALSE, fileEncoding = "UTF-8")

  plot_single(
    win_df = obj$win_df,
    summary_row = obj$summary,
    accession = acc,
    pretty_label = lab,
    outprefix = paste0(acc, "_window_genomic"),
    outdir = file.path(out_acc, "figures")
  )

  all_summaries[[acc]] <- obj$summary %>% mutate(accession = acc, label = lab)
  all_windows[[acc]] <- obj$win_df %>% mutate(accession = acc, label = lab)
}

SUM <- bind_rows(all_summaries)
WIN <- bind_rows(all_windows)

write.csv(SUM, file.path(OUTDIR, "ALL_window_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(WIN, file.path(OUTDIR, "ALL_window_CDS_table.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## ---------------------------
## 4) 多株相对坐标对齐总图（改为 facet，避免 y 轴拥挤）
## ---------------------------
labs_order <- manifest_ok$label
labs_order <- labs_order[!duplicated(labs_order)]

WIN_plot <- WIN %>%
  mutate(label = factor(label, levels = labs_order))

poly_list <- lapply(split(WIN_plot, WIN_plot$label), function(df) {
  df <- df %>% arrange(xmin_rel)
  if (nrow(df) == 0) return(NULL)
  make_arrow_polys(df, "xmin_rel", "xmax_rel", 1, ARROW_HEIGHT_Y, ARROW_HEAD_LEN_BP) %>%
    mutate(label = unique(df$label)[1])
})
poly_list <- poly_list[!vapply(poly_list, is.null, logical(1))]
POLY <- bind_rows(poly_list)

LABD <- WIN_plot %>%
  filter(label_show != "") %>%
  mutate(
    x = gene_center_rel,
    y = 1,
    label = factor(label, levels = labs_order)
  )

PLOT_HEIGHT_MULTI <- max(4.8, 1.05 * length(labs_order) + 0.9)

p_all <- ggplot(POLY, aes(x = x, y = y, group = gene_id, fill = fill_group)) +
  geom_polygon(color = "white", linewidth = 0.22) +
  geom_polygon(data = POLY %>% filter(is_focal), fill = NA, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = WINDOW_BP + 1, linetype = "dashed", color = "grey55", linewidth = 0.3) +
  geom_text_repel(
    data = LABD,
    aes(x = x, y = y, label = label_show),
    inherit.aes = FALSE,
    direction = "x",
    nudge_y = 0.28,
    box.padding = 0.22,
    point.padding = 0.08,
    segment.color = "grey55",
    segment.size = 0.25,
    min.segment.length = 0,
    max.overlaps = Inf,
    force = 0.8,
    size = LABEL_TEXT_SIZE,
    fontface = LABEL_FONTFACE,
    seed = 1
  ) +
  coord_cartesian(xlim = c(1, 2 * WINDOW_BP + 1), ylim = c(0.58, 1.72), clip = "off") +
  facet_grid(label ~ ., switch = "y") +
  labs(
    x = paste0("Relative position within window (bp) [", 2 * WINDOW_BP, " bp window]"),
    y = NULL,
    title = "Putative tdc/tyrDC regions (aligned by relative coordinates)"
  ) +
  base_theme_multi()

ggsave(file.path(OUTDIR, "ALL_tdc_operon_aligned_relative.pdf"), p_all,
       width = PLOT_WIDTH_MULTI, height = PLOT_HEIGHT_MULTI,
       device = "pdf", useDingbats = FALSE)
if (isTRUE(MAKE_PNG)) {
  ggsave(file.path(OUTDIR, "ALL_tdc_operon_aligned_relative.png"), p_all,
         width = PLOT_WIDTH_MULTI, height = PLOT_HEIGHT_MULTI, dpi = DPI_PNG)
}

sink(file.path(OUTDIR, "sessionInfo.txt"))
print(sessionInfo())
sink()

message("Done. Outputs saved in: ", normalizePath(OUTDIR))
