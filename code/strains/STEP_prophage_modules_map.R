local({

## =========================================================
## STEP_prophage_modules_map_FIXED_v5.R
## Prophage module annotation + module maps (robust, no unit errors)
##
## Inputs (in working dir or give full paths):
##   - prophage_sequence_stats.csv
##   - Efa13_chr.gff3, Efi9_chr.gff3
##   - Efa13_chr_swissprot.xls, Efi9_chr_swissprot.xls   (TSV text; may have # comment lines)
## Optional:
##   - crossStrain_nearestPairs_kmer4.csv
##
## Outputs (new folder):
##   - 01_prophage_genes_with_modules.csv
##   - 02_prophage_module_presence_summary.csv
##   - 03_best_hit_swissprot_per_orf.csv
##   - FigS_modules_composition_bar.pdf (+ optional PNG)
##   - FigS_modules_Efa13_prophage_maps.pdf
##   - FigS_modules_Efi9_prophage_maps.pdf
##   - FigS_modules_nearestPairs_topN.pdf (optional)
## =========================================================

## ---------------------------
## 0) Parameters (你主要改这里)
## ---------------------------
PROPH_STATS <- "prophage_sequence_stats.csv"

GFF_EFA <- "Efa13_chr.gff3"
GFF_EFI <- "Efi9_chr.gff3"

SWISS_EFA <- "Efa13_chr_swissprot.xls"
SWISS_EFI <- "Efi9_chr_swissprot.xls"

PAIR_FILE <- "crossStrain_nearestPairs_kmer4.csv"  # optional
TOP_N_PAIRS_TO_PLOT <- 5

## --- 标签策略：TRUE=只标关键模块基因(或GeneName)，更清爽；FALSE=尽量都标 ---
LABEL_KEY_GENES_ONLY <- TRUE

## --- 箭头/字体大小（核心）---
ARROW_HEIGHT <- 0.80       # 箭头“高度/粗细”（越大越高）
ARROW_HEAD_LEN_BP <- 280   # 箭头头在 x 轴长度（单位：bp，越大头越长）
HEAD_WIDTH_RATIO <- 1.40   # 箭头头部相对身体更“尖/更宽”（>1 更宽）

LABEL_SIZE <- 2          # 基因标签字号（geom_text 的 size）
TITLE_SIZE <- 12
AXIS_SIZE  <- 10
TRACK_LABEL_SIZE <- 10

## --- 图尺寸（inch）---
PDF_W_SINGLE <- 12
PDF_H_SINGLE <- 3.2
PDF_W_PAIR   <- 12
PDF_H_PAIR   <- 4.2

## 输出
OUTDIR <- paste0("OUT_prophage_modules_FIXEDv5_", format(Sys.time(), "%Y%m%d_%H%M%S"))
MAKE_PNG <- TRUE
DPI_PNG  <- 300

## ---------------------------
## 1) Packages
## ---------------------------
pkgs <- c("dplyr","ggplot2","readr","stringr","tibble","tidyr")
miss <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if(length(miss)>0) install.packages(miss, repos="https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "PNG"), showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(PROPH_STATS), file.exists(GFF_EFA), file.exists(GFF_EFI),
          file.exists(SWISS_EFA), file.exists(SWISS_EFI))

message("Output folder: ", normalizePath(OUTDIR))

## ---------------------------
## 2) Helpers
## ---------------------------

## 2.1 Read SwissProt TSV-with-#comments robustly (base R)
read_swiss_raw <- function(path){
  try_enc <- c("UTF-8-BOM","UTF-8","GB18030")
  last_err <- NULL
  for(enc in try_enc){
    x <- try(
      read.delim(
        file = path,
        sep = "\t",
        header = TRUE,
        comment.char = "#",
        quote = "",
        stringsAsFactors = FALSE,
        check.names = FALSE,
        fileEncoding = enc
      ),
      silent = TRUE
    )
    if(!inherits(x, "try-error")){
      x <- as_tibble(x)
      if(ncol(x) >= 5) return(x)
      last_err <- paste0("Read ok but too few columns (", ncol(x), ") with encoding=", enc)
    } else {
      last_err <- as.character(x)
    }
  }
  stop("Cannot read SwissProt TSV: ", path, "\nLast: ", last_err, call.=FALSE)
}

## 2.2 Normalize column names (handle full-width letters if stringi exists)
normalize_colnames <- function(nms){
  n0 <- nms
  if(requireNamespace("stringi", quietly = TRUE)){
    n0 <- stringi::stri_trans_nfkc(n0)
  }
  n1 <- tolower(gsub("[^a-z0-9]+", "_", n0))
  n1 <- gsub("^_+|_+$", "", n1)
  n1 <- gsub("_+", "_", n1)
  n1
}

find_col_regex <- function(nms, patterns){
  for(p in patterns){
    hit <- grep(p, nms, ignore.case = TRUE, value = TRUE)
    if(length(hit) > 0) return(hit[1])
  }
  NULL
}

## 2.3 Parse GFF3 CDS (no rtracklayer)
read_gff_cds <- function(gff, strain){
  g <- readr::read_tsv(
    gff,
    comment = "#",
    col_names = c("seqid","source","type","start","end","score","strand","phase","attributes"),
    col_types = readr::cols(
      seqid      = readr::col_character(),
      source     = readr::col_character(),
      type       = readr::col_character(),
      start      = readr::col_integer(),
      end        = readr::col_integer(),
      score      = readr::col_character(),
      strand     = readr::col_character(),
      phase      = readr::col_character(),
      attributes = readr::col_character()
    ),
    progress = FALSE
  )

  g <- g %>%
    mutate(type_lower = tolower(type)) %>%
    filter(type_lower == "cds") %>%
    mutate(
      strand = ifelse(strand %in% c("+","-"), strand, "+"),
      gene_len = end - start + 1,
      strain = strain
    ) %>%
    arrange(seqid, start, end)

  if(nrow(g) == 0) stop("No CDS entries found in GFF3: ", gff, call.=FALSE)

  parse_attr <- function(x, key){
    m <- str_match(x, paste0("(?i)(?:^|;)", key, "=([^;]+)"))
    out <- m[,2]
    ifelse(is.na(out), NA_character_, out)
  }

  g <- g %>%
    group_by(seqid) %>%
    mutate(
      orf_index = row_number(),
      ## 你的 SwissProt ORF 命名已验证：seqid_index (如 chr_1, chr_2 ...)
      orf_seqid_index = paste0(seqid, "_", orf_index),
      orf_attr_id      = parse_attr(attributes, "ID"),
      orf_attr_locus   = parse_attr(attributes, "locus_tag"),
      orf_attr_name    = parse_attr(attributes, "Name")
    ) %>%
    ungroup()

  g
}

## 2.4 SwissProt best-hit per ORF (supports name/escription/ene_ame)
read_swiss_besthit <- function(path){
  s <- read_swiss_raw(path)
  colnames(s) <- normalize_colnames(colnames(s))

  orf_col <- find_col_regex(colnames(s), c("^orf_name$","^orfname$","^qseqid$","^query$","^name$"))
  desc_col <- find_col_regex(colnames(s), c("^description$","description","^product$","product",
                                           "^annot$","annotation","^descr","descr",
                                           "^escription$","escription"))
  gene_col <- find_col_regex(colnames(s), c("^genename$","^gene_name$","^gene$","gene",
                                           "^ene_ame$","ene_ame"))
  evalue_col <- find_col_regex(colnames(s), c("^evalue$","^e_value$"))
  bits_col   <- find_col_regex(colnames(s), c("^bitscore$","^bit_score$","^score$"))
  swiss_id_col <- find_col_regex(colnames(s), c("^swissprot_id$","^sp_id$","^accession$","^sseqid$","swiss"))

  if(is.null(orf_col)) stop("Cannot find ORF name column in SwissProt file: ", path, call.=FALSE)
  if(is.null(evalue_col)) stop("SwissProt file missing evalue column: ", path, call.=FALSE)

  message("[SwissProt] ", basename(path),
          " | ORF=", orf_col,
          " | DESC=", ifelse(is.null(desc_col), "NA", desc_col),
          " | GENE=", ifelse(is.null(gene_col), "NA", gene_col))

  s2 <- s %>%
    transmute(
      orf_name = as.character(.data[[orf_col]]),
      evalue   = suppressWarnings(as.numeric(.data[[evalue_col]])),
      bitscore = if(!is.null(bits_col)) suppressWarnings(as.numeric(.data[[bits_col]])) else NA_real_,
      swissprot_id = if(!is.null(swiss_id_col)) as.character(.data[[swiss_id_col]]) else NA_character_,
      Description  = if(!is.null(desc_col)) as.character(.data[[desc_col]]) else NA_character_,
      GeneName     = if(!is.null(gene_col)) as.character(.data[[gene_col]]) else NA_character_
    ) %>%
    filter(!is.na(orf_name), nzchar(orf_name)) %>%
    group_by(orf_name) %>%
    arrange(evalue, desc(bitscore)) %>%
    slice(1) %>%
    ungroup()

  s2
}

choose_best_orfname <- function(gff_tbl, swiss_tbl){
  candidates <- c("orf_seqid_index","orf_attr_id","orf_attr_locus","orf_attr_name")
  swiss_set <- unique(swiss_tbl$orf_name)

  score_one <- function(col){
    v <- gff_tbl[[col]]
    sum(!is.na(v) & v %in% swiss_set)
  }
  scores <- sapply(candidates, score_one)
  best <- names(which.max(scores))

  message("ORF-name matching hits: ",
          paste0(names(scores), "=", scores, collapse=" | "),
          "  --> choose: ", best)
  best
}

## 2.5 Module classifier (keyword based)
classify_module <- function(desc, gene){
  txt <- paste(desc, gene, sep=" | ")
  txt <- tolower(ifelse(is.na(txt), "", txt))

  if(str_detect(txt, "integrase|site-specific recombinase|recombinase|tyrosine recombinase|serine recombinase|\\bint\\b")) return("Integrase")
  if(str_detect(txt, "repressor|\\bci\\b|\\bcro\\b|antirepressor|excisionase|\\bxis\\b|lysogen|regulator|transcriptional regulator")) return("Lysogeny/Regulation")
  if(str_detect(txt, "terminase|portal|capsid|head|prohead|scaffold|major capsid|minor capsid|packaging")) return("Packaging/Head")
  if(str_detect(txt, "tail|baseplate|tape measure|sheath|tube|tail fiber|tail spike")) return("Tail")
  if(str_detect(txt, "holin|endolysin|lysin|amidase|muramidase|lysozyme|peptidoglycan hydrolase|spanin")) return("Lysis")
  if(str_detect(txt, "helicase|primase|polymerase|\\bdna\\b|\\bssb\\b|replication|nuclease|exonuclease|ribonucleotide reductase")) return("Replication/Metabolism")

  if(is.na(desc) || desc=="" || str_detect(tolower(desc), "hypothetical|uncharacterized")) return("Other/Hypothetical")
  "Other"
}

module_short <- function(module){
  d <- c("Integrase"="int","Lysogeny/Regulation"="lys","Packaging/Head"="head","Tail"="tail",
         "Lysis"="lysis","Replication/Metabolism"="rep","Other/Hypothetical"="","Other"="")
  unname(d[module])
}

## 2.6 Intersect genes with prophage intervals
get_genes_in_prophage <- function(proph_tbl, gene_tbl){
  out <- vector("list", nrow(proph_tbl))
  for(i in seq_len(nrow(proph_tbl))){
    p <- proph_tbl[i,]
    g <- gene_tbl %>%
      filter(seqid == p$seq_id2, start <= p$stop, end >= p$start) %>%
      mutate(
        prophage_id = p$prophage_id,
        prophage_tag = p$prophage_tag,
        proph_start = p$start,
        proph_stop  = p$stop,
        rel_start = start - p$start,
        rel_end   = end   - p$start
      )
    out[[i]] <- g
  }
  bind_rows(out)
}

## 2.7 Manual arrow polygons (NO unit objects -> avoids upgradeUnit errors)
make_gene_polygons <- function(df, y_center=1, height=0.8, head_len_bp=250, head_width_ratio=1.4){
  if(nrow(df) == 0) return(tibble())

  df2 <- df %>%
    mutate(
      forward = strand == "+",
      len_bp = pmax(1, rel_end - rel_start),
      head_len = pmin(head_len_bp, len_bp * 0.45),
      body_y0 = y_center - height/2,
      body_y1 = y_center + height/2,
      head_y0 = y_center - (height/2) * head_width_ratio,
      head_y1 = y_center + (height/2) * head_width_ratio
    )

  poly_list <- lapply(seq_len(nrow(df2)), function(i){
    r <- df2[i,]
    if(isTRUE(r$forward)){
      x0 <- r$rel_start
      x2 <- r$rel_end
      x1 <- r$rel_end - r$head_len

      tibble(
        orf_name = r$orf_name,
        module = r$module,
        track = if("track" %in% names(r)) as.character(r$track) else NA_character_,
        y = c(r$body_y0, r$body_y0, r$head_y0, y_center, r$head_y1, r$body_y1, r$body_y1),
        x = c(x0, x1, x1, x2, x1, x1, x0),
        rel_start = r$rel_start,
        rel_end   = r$rel_end,
        strand    = r$strand,
        GeneName  = r$GeneName,
        Description = r$Description
      )
    } else {
      x0 <- r$rel_start
      x2 <- r$rel_end
      x1 <- r$rel_start + r$head_len

      tibble(
        orf_name = r$orf_name,
        module = r$module,
        track = if("track" %in% names(r)) as.character(r$track) else NA_character_,
        y = c(r$body_y0, r$body_y0, r$head_y0, y_center, r$head_y1, r$body_y1, r$body_y1),
        x = c(x2, x1, x1, x0, x1, x1, x2),
        rel_start = r$rel_start,
        rel_end   = r$rel_end,
        strand    = r$strand,
        GeneName  = r$GeneName,
        Description = r$Description
      )
    }
  })

  bind_rows(poly_list)
}

make_label_df <- function(df, y_center=1, key_only=TRUE){
  if(nrow(df) == 0) return(tibble())
  df %>%
    mutate(
      xmid = (rel_start + rel_end)/2,
      label_show = case_when(
        !is.na(GeneName) & GeneName != "" ~ GeneName,
        key_only & module %in% c("Integrase","Lysogeny/Regulation","Packaging/Head","Tail","Lysis") ~ module_short(module),
        TRUE ~ ""
      ),
      y = y_center
    ) %>%
    filter(label_show != "")
}

plot_one_prophage_manual <- function(df_one, pinfo){
  df_one <- df_one %>% arrange(rel_start)

  poly <- make_gene_polygons(df_one, y_center=1,
                            height=ARROW_HEIGHT,
                            head_len_bp=ARROW_HEAD_LEN_BP,
                            head_width_ratio=HEAD_WIDTH_RATIO)

  lab <- make_label_df(df_one, y_center=1, key_only=LABEL_KEY_GENES_ONLY)

  title_txt <- paste0(pinfo$strain, " ", pinfo$prophage_id,
                      " | ", round(pinfo$len_kb,2), " kb | GC ", round(pinfo$GC_pct,2), "% | ",
                      pinfo$seq_id2, ":", pinfo$start, "-", pinfo$stop)

  ggplot() +
    geom_polygon(data=poly, aes(x=x/1000, y=y, group=orf_name, fill=module), color="black", linewidth=0.25) +
    geom_text(data=lab, aes(x=xmid/1000, y=y, label=label_show), size=LABEL_SIZE, vjust=-0.8) +
    scale_y_continuous(limits=c(0.1, 1.9), breaks=NULL) +
    theme_bw() +
    theme(
      legend.position="bottom",
      plot.title = element_text(size=TITLE_SIZE, face="bold"),
      axis.title = element_text(size=AXIS_SIZE),
      axis.text  = element_text(size=AXIS_SIZE),
      panel.grid = element_blank()
    ) +
    labs(title=title_txt, x="Position within prophage (kb)", y=NULL, fill="Module")
}

plot_pair_manual <- function(df_pair, title_txt){
  # tracks: Efa13=2, Efi9=1
  df_pair <- df_pair %>% mutate(track = as.character(track))

  df_efa <- df_pair %>% filter(track=="Efa13") %>% arrange(rel_start)
  df_efi <- df_pair %>% filter(track=="Efi9") %>% arrange(rel_start)

  poly_efa <- make_gene_polygons(df_efa, y_center=2,
                                height=ARROW_HEIGHT,
                                head_len_bp=ARROW_HEAD_LEN_BP,
                                head_width_ratio=HEAD_WIDTH_RATIO)
  poly_efi <- make_gene_polygons(df_efi, y_center=1,
                                height=ARROW_HEIGHT,
                                head_len_bp=ARROW_HEAD_LEN_BP,
                                head_width_ratio=HEAD_WIDTH_RATIO)

  poly <- bind_rows(poly_efa, poly_efi)

  lab <- bind_rows(
    make_label_df(df_efa, y_center=2, key_only=TRUE),
    make_label_df(df_efi, y_center=1, key_only=TRUE)
  )

  xmax_kb <- max(df_pair$rel_end, na.rm=TRUE)/1000

  ggplot() +
    geom_polygon(data=poly, aes(x=x/1000, y=y, group=paste0(track,"__",orf_name), fill=module),
                 color="black", linewidth=0.25) +
    geom_text(data=lab, aes(x=xmid/1000, y=y, label=label_show), size=LABEL_SIZE, vjust=-0.8) +
    scale_y_continuous(
      breaks=c(1,2),
      labels=c("Efi9","Efa13"),
      limits=c(0.2, 2.8)
    ) +
    coord_cartesian(xlim=c(0, xmax_kb)) +
    theme_bw() +
    theme(
      legend.position="bottom",
      plot.title = element_text(size=TITLE_SIZE, face="bold"),
      axis.title = element_text(size=AXIS_SIZE),
      axis.text  = element_text(size=AXIS_SIZE),
      axis.text.y = element_text(size=TRACK_LABEL_SIZE, face="bold"),
      panel.grid = element_blank()
    ) +
    labs(title=title_txt, x="Position within prophage (kb)", y=NULL, fill="Module")
}

## ---------------------------
## 3) Load prophage regions table
## ---------------------------
proph <- readr::read_csv(PROPH_STATS, show_col_types = FALSE)
need_cols <- c("strain","prophage_id","prophage_tag","seq_id2","start","stop","len_kb","GC_pct")
miss_cols <- setdiff(need_cols, colnames(proph))
if(length(miss_cols)>0) stop("prophage_sequence_stats.csv missing: ", paste(miss_cols, collapse=", "), call.=FALSE)

proph <- proph %>%
  mutate(
    strain = as.character(strain),
    prophage_id = as.character(prophage_id),
    prophage_tag = as.character(prophage_tag),
    seq_id2 = as.character(seq_id2),
    start = as.integer(start),
    stop  = as.integer(stop)
  )

proph_efa <- proph %>% filter(strain=="Efa13")
proph_efi <- proph %>% filter(strain=="Efi9")

## ---------------------------
## 4) Load genes + SwissProt; join; classify modules
## ---------------------------
genes_efa <- read_gff_cds(GFF_EFA, "Efa13")
genes_efi <- read_gff_cds(GFF_EFI, "Efi9")

sp_efa <- read_swiss_besthit(SWISS_EFA)
sp_efi <- read_swiss_besthit(SWISS_EFI)

best_efa <- choose_best_orfname(genes_efa, sp_efa)
best_efi <- choose_best_orfname(genes_efi, sp_efi)

genes_efa <- genes_efa %>% mutate(orf_name = .data[[best_efa]])
genes_efi <- genes_efi %>% mutate(orf_name = .data[[best_efi]])

genes_efa2 <- genes_efa %>%
  left_join(sp_efa, by="orf_name") %>%
  mutate(
    Description = ifelse(is.na(Description) | Description=="", "hypothetical protein", Description),
    GeneName = ifelse(is.na(GeneName), "", GeneName),
    module = vapply(seq_along(Description), function(i) classify_module(Description[i], GeneName[i]), character(1))
  )

genes_efi2 <- genes_efi %>%
  left_join(sp_efi, by="orf_name") %>%
  mutate(
    Description = ifelse(is.na(Description) | Description=="", "hypothetical protein", Description),
    GeneName = ifelse(is.na(GeneName), "", GeneName),
    module = vapply(seq_along(Description), function(i) classify_module(Description[i], GeneName[i]), character(1))
  )

bind_rows(
  genes_efa2 %>% select(strain, orf_name, seqid, start, end, strand, swissprot_id, evalue, bitscore, GeneName, Description),
  genes_efi2 %>% select(strain, orf_name, seqid, start, end, strand, swissprot_id, evalue, bitscore, GeneName, Description)
) %>%
  write.csv(file.path(OUTDIR, "03_best_hit_swissprot_per_orf.csv"),
            row.names=FALSE, fileEncoding="UTF-8")

## ---------------------------
## 5) Intersect genes with prophage intervals
## ---------------------------
in_efa <- get_genes_in_prophage(proph_efa, genes_efa2) %>% mutate(track="Efa13")
in_efi <- get_genes_in_prophage(proph_efi, genes_efi2) %>% mutate(track="Efi9")

genes_in_prophage <- bind_rows(in_efa, in_efi) %>%
  mutate(forward = strand == "+") %>%
  select(strain, track, prophage_id, prophage_tag, seqid, proph_start, proph_stop,
         orf_name, start, end, strand, forward, rel_start, rel_end,
         module, GeneName, Description, swissprot_id, evalue, bitscore)

write.csv(genes_in_prophage, file.path(OUTDIR, "01_prophage_genes_with_modules.csv"),
          row.names=FALSE, fileEncoding="UTF-8")

## ---------------------------
## 6) Module presence summary + composition barplot
## ---------------------------
mod_sum <- genes_in_prophage %>%
  group_by(strain, prophage_tag, prophage_id, module) %>%
  summarise(n_genes = n(), .groups="drop") %>%
  tidyr::complete(module, nesting(strain, prophage_tag, prophage_id), fill=list(n_genes=0))

presence <- mod_sum %>%
  group_by(strain, prophage_tag, prophage_id) %>%
  summarise(
    Integrase = sum(n_genes[module=="Integrase"])>0,
    Lysogeny_Reg = sum(n_genes[module=="Lysogeny/Regulation"])>0,
    Head_Packaging = sum(n_genes[module=="Packaging/Head"])>0,
    Tail = sum(n_genes[module=="Tail"])>0,
    Lysis = sum(n_genes[module=="Lysis"])>0,
    Replication = sum(n_genes[module=="Replication/Metabolism"])>0,
    Other = sum(n_genes[module %in% c("Other","Other/Hypothetical")])>0,
    total_genes = sum(n_genes),
    .groups="drop"
  )

write.csv(presence, file.path(OUTDIR, "02_prophage_module_presence_summary.csv"),
          row.names=FALSE, fileEncoding="UTF-8")

p_bar <- ggplot(mod_sum, aes(x=prophage_tag, y=n_genes, fill=module)) +
  geom_col() +
  facet_wrap(~strain, scales="free_x", ncol=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        legend.position="bottom") +
  labs(x=NULL, y="Number of ORFs", title="Prophage module composition", fill="Module")

ggsave(file.path(OUTDIR, "FigS_modules_composition_bar.pdf"),
       p_bar, width=12, height=6, device="pdf", useDingbats=FALSE)
if(MAKE_PNG){
  ggsave(file.path(OUTDIR, "PNG", "FigS_modules_composition_bar.png"),
         p_bar, width=12, height=6, dpi=DPI_PNG)
}

## ---------------------------
## 7) Prophage module maps (one prophage per page) - MANUAL arrows
## ---------------------------
plot_prophage_set_manual <- function(strain_name, proph_tbl, genes_tbl, out_pdf){
  pdf(out_pdf, width=PDF_W_SINGLE, height=PDF_H_SINGLE, onefile=TRUE)
  on.exit(dev.off(), add=TRUE)

  for(i in seq_len(nrow(proph_tbl))){
    pinfo <- proph_tbl[i,]
    df_one <- genes_tbl %>% filter(prophage_tag == pinfo$prophage_tag) %>% arrange(rel_start)

    if(nrow(df_one)==0){
      plot.new()
      title(main=paste0(strain_name," ",pinfo$prophage_id," (no CDS overlapped)"))
      next
    }

    print(plot_one_prophage_manual(df_one, pinfo))
  }
}

plot_prophage_set_manual(
  "Efa13", proph_efa,
  genes_in_prophage %>% filter(strain=="Efa13"),
  file.path(OUTDIR, "FigS_modules_Efa13_prophage_maps.pdf")
)

plot_prophage_set_manual(
  "Efi9", proph_efi,
  genes_in_prophage %>% filter(strain=="Efi9"),
  file.path(OUTDIR, "FigS_modules_Efi9_prophage_maps.pdf")
)

## ---------------------------
## 8) Nearest-pair module comparison (optional)
## ---------------------------
if(file.exists(PAIR_FILE)){
  pairs <- readr::read_csv(PAIR_FILE, show_col_types = FALSE) %>%
    arrange(kmer_dist) %>%
    slice(1:TOP_N_PAIRS_TO_PLOT)

  pdf(file.path(OUTDIR, paste0("FigS_modules_nearestPairs_top",TOP_N_PAIRS_TO_PLOT,".pdf")),
      width=PDF_W_PAIR, height=PDF_H_PAIR, onefile=TRUE)

  for(i in seq_len(nrow(pairs))){
    p <- pairs[i,]

    efa_tag <- as.character(p$Efa13_tag)
    efi_tag <- as.character(p$Efi9_tag)

    df_efa <- genes_in_prophage %>% filter(prophage_tag == efa_tag) %>% mutate(track="Efa13")
    df_efi <- genes_in_prophage %>% filter(prophage_tag == efi_tag) %>% mutate(track="Efi9")
    df_pair <- bind_rows(df_efa, df_efi)

    if(nrow(df_pair)==0){
      plot.new(); title(main=paste0(efa_tag," vs ",efi_tag," (no genes overlapped)"))
      next
    }

    title_txt <- paste0(efa_tag," vs ",efi_tag," | kmer_dist=", signif(p$kmer_dist, 3))
    print(plot_pair_manual(df_pair, title_txt))
  }
  dev.off()
}

sink(file.path(OUTDIR, "sessionInfo.txt"))
print(sessionInfo())
sink()

message("Done. Outputs saved in: ", normalizePath(OUTDIR))

})
