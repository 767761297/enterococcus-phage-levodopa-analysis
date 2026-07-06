## 03_Composition_StackedBar_6Levels.R
## Stacked barplots for 6 taxonomic levels:
## phylum / class / order / family / genus / species
## Input: allsample_bundle_from_featuretable/allsample.<level>.reabundance.xls

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

# ------------------- SETTINGS -------------------
IN_DIR  <- "allsample_bundle_from_featuretable"
OUT_DIR <- "16S_Composition_StackedBars_6Levels"
dir.create(OUT_DIR, showWarnings = FALSE)

# sample order in your allsample tables
SAMPLE_ORDER <- c(paste0("A",1:5), paste0("C",1:5), paste0("D",1:5), paste0("B",1:5))

# group mapping: A=GR, B=MR, C=MR+Phage cocktail, D=MR+AFMT
map_group <- function(s) {
  if (str_detect(s, "^A")) return("GR")
  if (str_detect(s, "^B")) return("MR")
  if (str_detect(s, "^C")) return("MR+Phage cocktail")
  if (str_detect(s, "^D")) return("MR+AFMT")
  return(NA_character_)
}
GROUP_LEVELS <- c("GR","MR","MR+Phage cocktail","MR+AFMT")

# Levels to plot (must have corresponding allsample.<level>.reabundance.xls)
LEVELS <- c("phylum","class","order","family","genus","species")

# Keep only Top N taxa (others -> Other). Adjust as you like.
TOPN <- list(
  phylum  = 15,
  class   = 20,
  order   = 25,
  family  = 30,
  genus   = 25,
  species = 30
)

# Optional: minimum mean relative abundance to keep as named taxa
# (set to 0 to disable)
MIN_MEAN <- list(
  phylum  = 0,
  class   = 0,
  order   = 0,
  family  = 0,
  genus   = 0,
  species = 0
)

# Plot size
PLOT_W <- 12
PLOT_H <- 4.8
# ------------------------------------------------


# ------------------- HELPERS -------------------
read_reabundance <- function(level) {
  f <- file.path(IN_DIR, paste0("allsample.", level, ".reabundance.xls"))
  if (!file.exists(f)) stop("File not found: ", f)

  df <- read.delim(f, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(df) < 2) stop("Bad table (need >=2 columns): ", f)

  tax_col <- colnames(df)[1]

  df_long <- df %>%
    rename(taxon = all_of(tax_col)) %>%
    pivot_longer(cols = -taxon, names_to = "sample", values_to = "relab") %>%
    mutate(
      sample = as.character(sample),
      relab  = suppressWarnings(as.numeric(relab)),
      relab  = ifelse(is.na(relab), 0, relab),
      group  = vapply(sample, map_group, character(1))
    ) %>%
    filter(!is.na(group))

  # enforce order
  df_long$sample <- factor(df_long$sample, levels = SAMPLE_ORDER)
  df_long$group  <- factor(df_long$group, levels = GROUP_LEVELS)

  df_long
}

collapse_other <- function(df_long, topn, min_mean) {
  keep_tbl <- df_long %>%
    group_by(taxon) %>%
    summarise(mean_relab = mean(relab, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_relab)) %>%
    filter(mean_relab >= min_mean)

  keep_taxa <- head(keep_tbl$taxon, topn)

  df2 <- df_long %>%
    mutate(taxon2 = ifelse(taxon %in% keep_taxa, taxon, "Other")) %>%
    group_by(group, sample, taxon2) %>%
    summarise(relab = sum(relab, na.rm = TRUE), .groups = "drop") %>%
    rename(taxon = taxon2)

  # order legend by overall mean; Other last
  tax_order <- df2 %>%
    group_by(taxon) %>%
    summarise(m = mean(relab, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(m)) %>%
    pull(taxon)
  tax_order <- c(setdiff(tax_order, "Other"), "Other")
  df2$taxon <- factor(df2$taxon, levels = tax_order)

  df2
}

plot_stacked_samples <- function(level) {
  df <- read_reabundance(level)
  df <- collapse_other(df, topn = TOPN[[level]], min_mean = MIN_MEAN[[level]])

  p <- ggplot(df, aes(x = sample, y = relab, fill = taxon)) +
    geom_col(width = 0.95) +
    facet_grid(~ group, scales = "free_x", space = "free_x") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "Relative abundance", fill = NULL) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x        = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position    = "right",
      strip.background   = element_rect(fill = "white")
    )

  out_png <- file.path(OUT_DIR, paste0("StackedBar_", level, "_samples.png"))
  out_pdf <- file.path(OUT_DIR, paste0("StackedBar_", level, "_samples.pdf"))

  ggsave(out_png, p, width = PLOT_W, height = PLOT_H, dpi = 300)
  ggsave(out_pdf, p, width = PLOT_W, height = PLOT_H)

  message("[OK] saved: ", out_png)
  invisible(p)
}

# ------------------- RUN -------------------
# sanity check
missing <- LEVELS[!file.exists(file.path(IN_DIR, paste0("allsample.", LEVELS, ".reabundance.xls")))]
if (length(missing) > 0) {
  stop("Missing reabundance files for levels: ", paste(missing, collapse = ", "),
       "\nExpected under: ", IN_DIR)
}

for (lv in LEVELS) {
  plot_stacked_samples(lv)
}

message("[DONE] All 6 stacked barplots saved in: ", OUT_DIR)
