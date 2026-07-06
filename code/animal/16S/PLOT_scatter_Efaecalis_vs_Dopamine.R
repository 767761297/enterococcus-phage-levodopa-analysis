library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)

## ---- inputs ----
SPECIES_DIRS <- c(".", "allsample_bundle_from_featuretable")
SPECIES_FILE_PAT <- "allsample\\.species\\.(reabundance|abundance)\\.xls$"
MET_QUANT <- "Q2.xlsx"

TARGET_TAXON <- "Enterococcus_faecalis"
TARGET_MET   <- "Dopamine"

DROP_AFMT <- TRUE   # 若要剔除D组，把这里改 TRUE; 不剔除D组则为FALSE

## ---- helper: find species table ----
find_species_file <- function() {
  hits <- c()
  for (d in SPECIES_DIRS) {
    if (!dir.exists(d) && d != ".") next
    hits <- c(hits, list.files(d, pattern = SPECIES_FILE_PAT, recursive = TRUE,
                               full.names = TRUE, ignore.case = TRUE))
  }
  hits <- unique(hits)
  if (length(hits) == 0) stop("Cannot find allsample.species.reabundance/abundance table.")
  # prefer reabundance
  reab <- hits[grepl("reabundance", hits, ignore.case = TRUE)]
  if (length(reab) > 0) return(reab[1])
  hits[1]
}

infer_group <- function(s) {
  dplyr::case_when(
    str_detect(s, "^A\\d+$") ~ "GR",
    str_detect(s, "^B\\d+$") ~ "MR",
    str_detect(s, "^C\\d+$") ~ "MR+Phage cocktail",
    str_detect(s, "^D\\d+$") ~ "MR+AFMT",
    TRUE ~ NA_character_
  )
}

clr_log2 <- function(mat) {
  x <- as.matrix(mat)
  x[is.na(x)] <- 0
  x[x < 0] <- 0
  for (j in seq_len(ncol(x))) {
    v <- x[, j]
    pos <- v[v > 0]
    pc <- if (length(pos) > 0) max(min(pos)/2, 1e-12) else 1e-6
    v2 <- v + pc
    gm <- exp(mean(log(v2)))
    x[, j] <- log2(v2 / gm)
  }
  x
}

## ---- 1) read species abundance ----
sp_file <- find_species_file()
df_sp <- read.delim(sp_file, sep = "\t", header = TRUE, check.names = FALSE,
                    stringsAsFactors = FALSE, quote = "", comment.char = "")

sample_cols <- grep("^[ABCD][1-5]$", names(df_sp), value = TRUE)
if (length(sample_cols) != 20) stop("Cannot detect 20 sample columns A1..D5 in species table.")

if (DROP_AFMT) sample_cols <- sample_cols[!str_detect(sample_cols, "^D")]

tax_col <- names(df_sp)[1]
mat_sp <- df_sp %>%
  transmute(taxon = as.character(.data[[tax_col]]) %>% str_trim()) %>%
  bind_cols(df_sp[, sample_cols, drop = FALSE]) %>%
  filter(!is.na(taxon) & taxon != "")

mat <- as.matrix(mat_sp[, sample_cols, drop = FALSE])
mat <- apply(mat, 2, function(x) suppressWarnings(as.numeric(x)))
rownames(mat) <- mat_sp$taxon

# CLR
mat_clr <- clr_log2(mat)

if (!(TARGET_TAXON %in% rownames(mat_clr))) stop("Taxon not found in species table: ", TARGET_TAXON)
x <- mat_clr[TARGET_TAXON, ]

## ---- 2) read metabolomics quant (Dopamine) ----
met_raw <- read_excel(MET_QUANT, sheet = excel_sheets(MET_QUANT)[1]) %>% as.data.frame(check.names = FALSE)
met_smp <- grep("^[ABCD][1-5]$", names(met_raw), value = TRUE)
if (length(met_smp) != 20) stop("Cannot detect 20 sample columns A1..D5 in metabolomics quant.")

if (DROP_AFMT) met_smp <- met_smp[!str_detect(met_smp, "^D")]

compound_col <- grep("Compound\\s*Name", names(met_raw), ignore.case = TRUE, value = TRUE)[1]
if (is.na(compound_col)) stop("Cannot find Compound Name column in metabolomics quant.")

met_raw$metabolite <- as.character(met_raw[[compound_col]]) %>% str_trim()
row <- met_raw %>% filter(metabolite == TARGET_MET)
if (nrow(row) == 0) stop("Metabolite not found in quant: ", TARGET_MET)

y0 <- as.numeric(row[1, met_smp, drop = TRUE])
# log2 with pseudocount
pos <- y0[y0 > 0]
pc <- if (length(pos) > 0) max(min(pos)/2, 1e-12) else 1e-6
y <- log2(y0 + pc)

## ---- 3) build plot df ----
plot_df <- data.frame(
  sample = met_smp,
  group = factor(sapply(met_smp, infer_group),
                 levels = c("GR","MR","MR+Phage cocktail","MR+AFMT")),
  Enterococcus_faecalis_CLR = as.numeric(x[met_smp]),
  Dopamine_log2 = as.numeric(y),
  stringsAsFactors = FALSE
)

plot_df <- plot_df %>% filter(!is.na(group))
if (DROP_AFMT) plot_df <- plot_df %>% filter(group != "MR+AFMT")

## ---- 4) stats ----
ct <- suppressWarnings(cor.test(plot_df$Enterococcus_faecalis_CLR,
                               plot_df$Dopamine_log2,
                               method = "spearman", exact = FALSE))
rho <- unname(ct$estimate); p <- ct$p.value; n <- nrow(plot_df)

## ---- 5) figure ----
p1 <- ggplot(plot_df, aes(x = Enterococcus_faecalis_CLR, y = Dopamine_log2, color = group)) +
  geom_point(size = 2.5, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.6) +
  theme_bw(base_size = 12) +
  labs(x = "Enterococcus faecalis (CLR, species level)",
       y = "Dopamine (log2(nmol/g + pseudocount))") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.2,
           label = sprintf("Spearman rho=%.3f\np=%.2e\nn=%d", rho, p, n))

ggsave(filename = file.path("9. Corr_3groups_noAFMT_16Ssig_q0.1_vs_Metsig_q0.2", 
                            ifelse(DROP_AFMT, "Scatter_Efaecalis_Dopamine_noAFMT.pdf",
                                           "Scatter_Efaecalis_Dopamine_all4groups.pdf")),
       plot = p1, width = 6.5, height = 5.5)
