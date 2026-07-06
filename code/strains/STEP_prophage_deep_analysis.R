local({

## =========================================================
## STEP_prophage_deep_analysis_FIXED2.R
## Deep prophage analysis using:
##   - prophage coordinate tables (*.xls, tab-delimited text)
##   - genome FASTA for each strain
##
## Outputs (NEW folder):
##   - Extracted prophage FASTA (per strain + combined)
##   - CSV tables: regions, sequence stats, kmer distance, nearest pairs, summary
##   - PDF figures (required): genome map, length dist, GC vs length, kmer PCA, distance heatmap
## =========================================================

## ---------------------------
## 0) Parameters (你主要改这里)
## ---------------------------
EFA_PH <- "Efa13_chr_prophage.xls"
EFI_PH <- "Efi9_chr_prophage.xls"

EFA_FA <- "Efa13_.fasta"
EFI_FA <- "Efi9_chr.fasta"

KMER_K <- 4

MAKE_PNG <- TRUE
DPI_PNG  <- 300

OUTDIR <- paste0("OUT_prophage_deep_FIXED2_", format(Sys.time(), "%Y%m%d_%H%M%S"))

## ---------------------------
## 1) Packages
## ---------------------------
pkgs <- c("dplyr","ggplot2","readr","stringr","tibble")
miss <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if(length(miss)>0) install.packages(miss, repos="https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(readr); library(stringr); library(tibble)
})

## ---------------------------
## 2) Safety checks
## ---------------------------
must_string <- function(x, nm){
  if(!is.character(x) || length(x)!=1 || is.na(x) || !nzchar(x)){
    stop("Parameter ", nm, " must be a non-empty character string. class=",
         paste(class(x), collapse=","), call.=FALSE)
  }
  TRUE
}
must_string(EFA_PH,"EFA_PH"); must_string(EFI_PH,"EFI_PH")
must_string(EFA_FA,"EFA_FA"); must_string(EFI_FA,"EFI_FA")
must_string(OUTDIR,"OUTDIR")

stopifnot(file.exists(EFA_PH), file.exists(EFI_PH), file.exists(EFA_FA), file.exists(EFI_FA))

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "FASTA_prophage"), showWarnings = FALSE, recursive = TRUE)
message("Output folder: ", normalizePath(OUTDIR))

## ---------------------------
## 3) Helpers
## ---------------------------

read_fasta_named <- function(fa_file){
  x <- readLines(fa_file, warn = FALSE, encoding="UTF-8")
  x <- x[nzchar(x)]
  hdr <- grep("^>", x)
  if(length(hdr)==0) stop("No FASTA headers found: ", fa_file, call.=FALSE)

  ids <- sub("^>\\s*", "", x[hdr])
  ids <- sub("\\s.*$", "", ids)
  seqs <- character(length(ids))

  for(i in seq_along(hdr)){
    s <- hdr[i] + 1
    e <- if(i < length(hdr)) hdr[i+1]-1 else length(x)
    seqs[i] <- paste0(gsub("\\s+", "", x[s:e]), collapse="")
  }
  names(seqs) <- ids
  seqs
}

read_prophage_tbl <- function(file, strain){
  lines <- readLines(file, warn=FALSE, encoding="UTF-8")
  lines <- lines[!grepl("^#", lines)]
  lines <- lines[nzchar(lines)]
  spl <- strsplit(lines, "\t", fixed=TRUE)

  dat <- spl[-1]  # skip header
  if(length(dat)==0) stop("No data rows in: ", file, call.=FALSE)

  maxlen <- max(lengths(dat), 10L)
  base_names <- c("prophage_id","seq_id","start","stop",
                  "attL_start","attL_end","attR_start","attR_end",
                  "attL_seq","attR_seq")
  coln <- base_names[seq_len(maxlen)]

  dat2 <- lapply(dat, function(x){ length(x) <- maxlen; x })
  df <- as.data.frame(do.call(rbind, dat2), stringsAsFactors=FALSE)
  colnames(df) <- coln

  num_cols <- c("start","stop","attL_start","attL_end","attR_start","attR_end")
  for(nm in intersect(num_cols, colnames(df))){
    df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
  }

  df %>%
    mutate(
      strain = strain,
      seq_id = sub("\\s.*$", "", seq_id),
      len_bp = stop - start + 1,
      len_kb = len_bp/1000,
      attL_len = ifelse(!is.na(attL_seq) & nzchar(attL_seq), nchar(attL_seq), NA_integer_),
      attR_len = ifelse(!is.na(attR_seq) & nzchar(attR_seq), nchar(attR_seq), NA_integer_),
      att_same = ifelse(!is.na(attL_seq) & !is.na(attR_seq) & attL_seq == attR_seq, TRUE, FALSE)
    )
}

resolve_contig <- function(id, genome_names){
  if(id %in% genome_names) return(id)

  hit <- genome_names[str_detect(genome_names, fixed(id))]
  if(length(hit)==1) return(hit[1])

  strip <- function(z){
    z <- tolower(z)
    z <- gsub("^chr[_-]*", "", z)
    z <- gsub("^chromosome[_-]*", "", z)
    z <- gsub("[^a-z0-9]+", "", z)
    z
  }
  id3 <- strip(id); gn3 <- strip(genome_names)
  hit3 <- genome_names[gn3 == id3]
  if(length(hit3)==1) return(hit3[1])

  stop("Cannot match prophage seq_id '", id, "' to FASTA contig names.\nFASTA contigs: ",
       paste(genome_names, collapse=", "), call.=FALSE)
}

extract_seq <- function(genome_named, contig, start, stop){
  g <- genome_named[[contig]]
  L <- nchar(g)
  if(is.na(start) || is.na(stop) || start < 1 || stop < 1) return(NA_character_)
  s <- max(1, start); e <- min(L, stop)
  if(s > e) return(NA_character_)
  substr(g, s, e)
}

write_fasta <- function(headers, seqs, out_file){
  con <- file(out_file, open="wt", encoding="UTF-8")
  on.exit(close(con), add=TRUE)
  for(i in seq_along(seqs)){
    writeLines(paste0(">", headers[i]), con)
    ss <- seqs[i]
    if(is.na(ss) || !nzchar(ss)){ writeLines("N", con); next }
    ss <- gsub("\\s+", "", ss)
    chunks <- substring(ss, seq(1, nchar(ss), 80), pmin(nchar(ss), seq(1, nchar(ss), 80)+79))
    writeLines(chunks, con)
  }
}

## --- vectorized seq stats ---
calc_one_stats <- function(s){
  if(is.na(s) || !nzchar(s)) return(c(len_bp2=NA_real_, GC_pct=NA_real_, N_pct=NA_real_))
  x <- toupper(as.character(s))
  x <- gsub("[^ACGTN]", "N", x)
  n <- nchar(x)
  A <- str_count(x, "A"); C <- str_count(x, "C"); G <- str_count(x, "G"); T <- str_count(x, "T")
  Nn <- str_count(x, "N")
  denom <- max(1, (A+C+G+T))
  c(len_bp2=as.numeric(n), GC_pct=100*(G+C)/denom, N_pct=100*Nn/max(1,n))
}
calc_stats_vec <- function(seq_vec){
  mat <- t(vapply(seq_vec, calc_one_stats, numeric(3)))
  as_tibble(mat)
}

## ---------------------------
## 4) Load
## ---------------------------
efa_tbl <- read_prophage_tbl(EFA_PH, "Efa13")
efi_tbl <- read_prophage_tbl(EFI_PH, "Efi9")

efa_gen <- read_fasta_named(EFA_FA)
efi_gen <- read_fasta_named(EFI_FA)

efa_tbl <- efa_tbl %>% mutate(seq_id2 = vapply(seq_id, resolve_contig, character(1), genome_names=names(efa_gen)))
efi_tbl <- efi_tbl %>% mutate(seq_id2 = vapply(seq_id, resolve_contig, character(1), genome_names=names(efi_gen)))

all_tbl <- bind_rows(efa_tbl, efi_tbl) %>% arrange(strain, seq_id2, start)
write.csv(all_tbl, file.path(OUTDIR, "prophage_regions_all.csv"), row.names=FALSE, fileEncoding="UTF-8")

## ---------------------------
## 5) Extract sequences
## ---------------------------
efa_seqs <- mapply(function(ctg, s, e) extract_seq(efa_gen, ctg, s, e),
                   efa_tbl$seq_id2, efa_tbl$start, efa_tbl$stop, USE.NAMES=FALSE)
efi_seqs <- mapply(function(ctg, s, e) extract_seq(efi_gen, ctg, s, e),
                   efi_tbl$seq_id2, efi_tbl$start, efi_tbl$stop, USE.NAMES=FALSE)

efa_headers <- paste0("Efa13|", efa_tbl$prophage_id, "|", efa_tbl$seq_id2, ":", efa_tbl$start, "-", efa_tbl$stop)
efi_headers <- paste0("Efi9|",  efi_tbl$prophage_id, "|", efi_tbl$seq_id2, ":", efi_tbl$start, "-", efi_tbl$stop)

write_fasta(efa_headers, efa_seqs, file.path(OUTDIR, "FASTA_prophage", "Efa13_prophage_extracted.fasta"))
write_fasta(efi_headers, efi_seqs, file.path(OUTDIR, "FASTA_prophage", "Efi9_prophage_extracted.fasta"))
write_fasta(c(efa_headers, efi_headers), c(efa_seqs, efi_seqs),
            file.path(OUTDIR, "FASTA_prophage", "BOTH_prophage_extracted.fasta"))

## ---------------------------
## 6) Sequence stats (✅ FIX: do not reference seq_str inside pipe)
## ---------------------------
efa_seq_str <- as.character(efa_seqs)
efi_seq_str <- as.character(efi_seqs)

efa_stats <- bind_cols(
  efa_tbl,
  tibble(seq_str = efa_seq_str),
  calc_stats_vec(efa_seq_str)
) %>% mutate(prophage_tag = paste0(strain, "_", prophage_id))

efi_stats <- bind_cols(
  efi_tbl,
  tibble(seq_str = efi_seq_str),
  calc_stats_vec(efi_seq_str)
) %>% mutate(prophage_tag = paste0(strain, "_", prophage_id))

stats_all <- bind_rows(efa_stats, efi_stats) %>%
  mutate(
    len_bp = ifelse(!is.na(len_bp2), len_bp2, len_bp),
    len_kb = len_bp/1000
  ) %>%
  select(strain, prophage_id, prophage_tag, seq_id2, start, stop,
         len_bp, len_kb, GC_pct, N_pct,
         attL_start, attL_end, attR_start, attR_end, attL_len, attR_len, att_same, attL_seq, attR_seq)

write.csv(stats_all, file.path(OUTDIR, "prophage_sequence_stats.csv"),
          row.names=FALSE, fileEncoding="UTF-8")

## ---------------------------
## 7) k-mer analysis
## ---------------------------
make_all_kmers <- function(k){
  alph <- c("A","C","G","T"); km <- alph
  if(k==1) return(km)
  for(i in 2:k) km <- as.vector(outer(km, alph, paste0))
  km
}
all_kmers <- make_all_kmers(KMER_K)

kmer_freq <- function(s, k, levels_km){
  if(is.na(s) || !nzchar(s)) return(rep(0, length(levels_km)))
  x <- toupper(as.character(s))
  x <- gsub("[^ACGT]", "N", x)
  n <- nchar(x)
  if(n < k) return(rep(0, length(levels_km)))
  km <- substring(x, 1:(n-k+1), k:n)
  km <- km[!grepl("N", km, fixed=TRUE)]
  if(length(km)==0) return(rep(0, length(levels_km)))
  tab <- table(factor(km, levels=levels_km))
  as.numeric(tab) / sum(tab)
}

seq_vec <- c(efa_seqs, efi_seqs)
tags <- c(paste0("Efa13_", efa_tbl$prophage_id), paste0("Efi9_", efi_tbl$prophage_id))
strain_vec <- c(rep("Efa13", length(efa_seqs)), rep("Efi9", length(efi_seqs)))

kmat <- t(vapply(seq_vec, function(s) kmer_freq(s, KMER_K, all_kmers),
                 numeric(length(all_kmers))))
rownames(kmat) <- tags

d <- dist(kmat, method="euclidean")
dmat <- as.matrix(d)
write.csv(dmat, file.path(OUTDIR, paste0("kmer",KMER_K,"_distance_matrix.csv")),
          fileEncoding="UTF-8")

efa_idx <- which(strain_vec=="Efa13")
efi_idx <- which(strain_vec=="Efi9")

nearest_pairs <- bind_rows(lapply(efa_idx, function(i){
  di <- dmat[i, efi_idx, drop=TRUE]
  j <- efi_idx[which.min(di)]
  tibble(Efa13_tag = rownames(dmat)[i],
         Efi9_tag  = rownames(dmat)[j],
         kmer_dist = min(di))
})) %>% arrange(kmer_dist)

write.csv(nearest_pairs, file.path(OUTDIR, paste0("crossStrain_nearestPairs_kmer",KMER_K,".csv")),
          row.names=FALSE, fileEncoding="UTF-8")

## ---------------------------
## 8) Figures (PDF required)
## ---------------------------
map_df <- all_tbl %>% mutate(start_mb=start/1e6, stop_mb=stop/1e6)
p1 <- ggplot(map_df) +
  geom_rect(aes(xmin=start_mb, xmax=stop_mb, ymin=0, ymax=1),
            fill="grey20", alpha=0.75) +
  facet_wrap(~strain, ncol=1, scales="free_x") +
  theme_bw() +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
        strip.text=element_text(face="bold")) +
  labs(x="Genomic position (Mb)", y=NULL, title="Predicted prophage regions on genomes")
ggsave(file.path(OUTDIR, "Fig1_prophage_genomeMap_Mb.pdf"), p1, width=10, height=5,
       device="pdf", useDingbats=FALSE)
if(MAKE_PNG) ggsave(file.path(OUTDIR, "Fig1_prophage_genomeMap_Mb.png"), p1, width=10, height=5, dpi=DPI_PNG)

p2 <- ggplot(stats_all, aes(x=strain, y=len_kb)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.12, height=0, size=2, alpha=0.85) +
  theme_bw() +
  labs(x=NULL, y="Prophage length (kb)", title="Prophage length distribution")
ggsave(file.path(OUTDIR, "Fig2_prophage_length_distribution.pdf"), p2, width=6.5, height=4.5,
       device="pdf", useDingbats=FALSE)
if(MAKE_PNG) ggsave(file.path(OUTDIR, "Fig2_prophage_length_distribution.png"), p2, width=6.5, height=4.5, dpi=DPI_PNG)

p3 <- ggplot(stats_all, aes(x=len_kb, y=GC_pct, shape=strain)) +
  geom_point(size=3, alpha=0.9) +
  theme_bw() +
  labs(x="Prophage length (kb)", y="GC content (%)", title="Prophage length vs GC%")
ggsave(file.path(OUTDIR, "Fig3_length_vs_GC.pdf"), p3, width=7, height=5,
       device="pdf", useDingbats=FALSE)
if(MAKE_PNG) ggsave(file.path(OUTDIR, "Fig3_length_vs_GC.png"), p3, width=7, height=5, dpi=DPI_PNG)

pc <- prcomp(kmat, center=TRUE, scale.=FALSE)
pc_df <- tibble(tag=rownames(kmat), strain=strain_vec, PC1=pc$x[,1], PC2=pc$x[,2])
p4 <- ggplot(pc_df, aes(x=PC1, y=PC2, shape=strain)) +
  geom_point(size=3, alpha=0.9) +
  theme_bw() +
  labs(x="PC1 (k-mer profile)", y="PC2 (k-mer profile)",
       title=paste0("PCA of prophages by ", KMER_K, "-mer composition"))
ggsave(file.path(OUTDIR, paste0("Fig4_kmer",KMER_K,"_PCA.pdf")), p4, width=7, height=5,
       device="pdf", useDingbats=FALSE)
if(MAKE_PNG) ggsave(file.path(OUTDIR, paste0("Fig4_kmer",KMER_K,"_PCA.png")), p4, width=7, height=5, dpi=DPI_PNG)

hc <- hclust(d, method="average")
pdf(file.path(OUTDIR, paste0("Fig5_kmer",KMER_K,"_distance_heatmap.pdf")), width=9, height=9)
heatmap(dmat, Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc),
        scale="none", margins=c(7,7),
        main=paste0(KMER_K, "-mer Euclidean distance (prophage sequences)"))
dev.off()

## ---------------------------
## 9) Summary
## ---------------------------
summary_df <- stats_all %>%
  group_by(strain) %>%
  summarise(
    n_prophage = n(),
    total_kb = sum(len_kb, na.rm=TRUE),
    median_kb = median(len_kb, na.rm=TRUE),
    mean_GC = mean(GC_pct, na.rm=TRUE),
    mean_Npct = mean(N_pct, na.rm=TRUE),
    att_same_n = sum(att_same, na.rm=TRUE),
    .groups="drop"
  )
write.csv(summary_df, file.path(OUTDIR, "prophage_summary_by_strain.csv"),
          row.names=FALSE, fileEncoding="UTF-8")

sink(file.path(OUTDIR, "sessionInfo.txt"))
print(sessionInfo())
sink()

print(summary_df)
message("Done. All outputs saved in: ", normalizePath(OUTDIR))

})
