# STEP_make_allsample_like_bundle_from_featuretable.R
# Input (unchanged): featuretable_no_singletons.txt
# Output: "allsample.*.xls" tables in the same style as your example:
#   - allsample.feature.xls / allsample.feature.reabundance.xls
#   - allsample.taxonomy.xls (keeps the first comment line)
#   - allsample.all.abundance.xls / allsample.all.reabundance.xls
#   - allsample.full.abundance.xls / allsample.full.reabundance.xls
#   - allsample.{phylum,class,order,family,genus,species}.{abundance,reabundance}.xls
#
# Note: these ".xls" files are TAB-separated text (TSV), consistent with your current files.

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

# -------------------- SETTINGS --------------------
IN_FILE <- "featuretable_no_singletons.txt"
OUT_DIR <- "allsample_bundle_from_featuretable"

# Match your example column order: A -> C -> D -> B
SAMPLE_ORDER <- c(paste0("A", 1:5), paste0("C", 1:5), paste0("D", 1:5), paste0("B", 1:5))

REL_DIGITS <- 10            # printing digits for relative abundance
DROP_UNCLASSIFIED <- FALSE  # set TRUE if you want to remove "Unclassified" taxa rows
# --------------------------------------------------

dir.create(OUT_DIR, showWarnings = FALSE)
stopifnot(file.exists(IN_FILE))

# ---------- read featuretable (robust) ----------
skip_n <- 0
first_line <- readLines(IN_FILE, n = 1, warn = FALSE)
if (length(first_line) == 1 && grepl("^#\\s*Constructed", first_line)) skip_n <- 1

df <- read.delim(
  IN_FILE, sep = "\t", header = TRUE, skip = skip_n,
  check.names = FALSE, stringsAsFactors = FALSE,
  quote = "", comment.char = ""
)

colnames(df)[1] <- sub("^#", "", colnames(df)[1]) # "#FeatureID" -> "FeatureID"
stopifnot("FeatureID" %in% colnames(df))
stopifnot("taxonomy" %in% colnames(df))

# detect sample columns
sample_cols <- setdiff(colnames(df), c("FeatureID", "taxonomy"))
sample_cols <- sample_cols[grepl("^[A-D][0-9]+$", sample_cols)]
stopifnot(length(sample_cols) > 0)

# reorder samples (fallback to detected order if SAMPLE_ORDER not fully present)
if (all(SAMPLE_ORDER %in% sample_cols)) {
  sample_cols <- SAMPLE_ORDER
} else {
  message("[WARN] SAMPLE_ORDER not fully found; using detected sample order.")
}

# numeric counts
df[sample_cols] <- lapply(df[sample_cols], function(x) suppressWarnings(as.numeric(x)))
df[sample_cols] <- lapply(df[sample_cols], function(x) ifelse(is.na(x), 0, x))

# ---------- taxonomy parsing ----------
split_tax <- function(x) {
  parts <- unlist(strsplit(as.character(x), ";"))
  parts <- trimws(parts)
  parts[parts != ""]
}
get_rank <- function(tax, prefix) {
  parts <- split_tax(tax)
  hit <- parts[startsWith(parts, prefix)]
  if (length(hit) == 0) return(NA_character_)
  val <- sub(paste0("^", prefix), "", hit[1])
  val <- trimws(val)
  if (val == "" || tolower(val) == "na") return(NA_character_)
  val
}
ranks <- tibble(
  Kingdom = vapply(df$taxonomy, get_rank, character(1), prefix = "k__"),
  Phylum  = vapply(df$taxonomy, get_rank, character(1), prefix = "p__"),
  Class   = vapply(df$taxonomy, get_rank, character(1), prefix = "c__"),
  Order   = vapply(df$taxonomy, get_rank, character(1), prefix = "o__"),
  Family  = vapply(df$taxonomy, get_rank, character(1), prefix = "f__"),
  Genus   = vapply(df$taxonomy, get_rank, character(1), prefix = "g__"),
  Species = vapply(df$taxonomy, get_rank, character(1), prefix = "s__")
)
# fill missing
for (nm in names(ranks)) ranks[[nm]][is.na(ranks[[nm]])] <- "Unclassified"
# species fallback
idx <- ranks$Species == "Unclassified" & ranks$Genus != "Unclassified"
ranks$Species[idx] <- paste0(ranks$Genus[idx], "_sp")

# ---------- helpers ----------
write_tsv_xls <- function(dat, path) {
  write.table(dat, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

to_relab <- function(mat) {
  denom <- colSums(mat)
  denom[denom == 0] <- NA_real_
  rel <- sweep(mat, 2, denom, "/")
  rel[is.na(rel)] <- 0
  round(rel, REL_DIGITS)
}

agg_by <- function(key_cols, counts_df) {
  tmp <- bind_cols(ranks[key_cols], counts_df)
  if (DROP_UNCLASSIFIED) {
    keep <- apply(tmp[key_cols], 1, function(x) all(x != "Unclassified"))
    tmp <- tmp[keep, , drop = FALSE]
  }
  out <- tmp %>%
    group_by(across(all_of(key_cols))) %>%
    summarise(across(all_of(sample_cols), ~sum(.x, na.rm = TRUE)), .groups = "drop") %>%
    arrange(across(all_of(key_cols)))
  out
}

# ---------- 1) feature tables ----------
feature_ab <- df %>%
  select(FeatureID, all_of(sample_cols)) %>%
  rename(Feature = FeatureID)

write_tsv_xls(feature_ab, file.path(OUT_DIR, "allsample.feature.xls"))

feature_rea <- feature_ab
feature_rea[sample_cols] <- to_relab(as.matrix(feature_ab[sample_cols]))
write_tsv_xls(feature_rea, file.path(OUT_DIR, "allsample.feature.reabundance.xls"))

# ---------- 2) taxonomy table (with first comment line + #FeatureID header) ----------
tax_path <- file.path(OUT_DIR, "allsample.taxonomy.xls")
con <- file(tax_path, open = "wt", encoding = "UTF-8")
writeLines("# Constructed from biom file", con = con)

tax_tbl <- df %>% select(FeatureID, all_of(sample_cols), taxonomy)
colnames(tax_tbl)[1] <- "#FeatureID"
write_tsv_xls(tax_tbl, tax_path)  # overwrite without comment
# re-write with comment + body
con <- file(tax_path, open = "wt", encoding = "UTF-8")
writeLines("# Constructed from biom file", con = con)
close(con)
suppressWarnings(write.table(tax_tbl, tax_path, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE, append=TRUE))

# ---------- 3) all/full tables ----------
all_cols  <- c("Phylum","Class","Order","Family","Genus","Species")
full_cols <- c("Kingdom", all_cols)

all_ab <- agg_by(all_cols, df %>% select(all_of(sample_cols)))
write_tsv_xls(all_ab, file.path(OUT_DIR, "allsample.all.abundance.xls"))
all_rea <- all_ab
all_rea[sample_cols] <- to_relab(as.matrix(all_ab[sample_cols]))
write_tsv_xls(all_rea, file.path(OUT_DIR, "allsample.all.reabundance.xls"))

full_ab <- agg_by(full_cols, df %>% select(all_of(sample_cols)))
write_tsv_xls(full_ab, file.path(OUT_DIR, "allsample.full.abundance.xls"))
full_rea <- full_ab
full_rea[sample_cols] <- to_relab(as.matrix(full_ab[sample_cols]))
write_tsv_xls(full_rea, file.path(OUT_DIR, "allsample.full.reabundance.xls"))

# ---------- 4) rank-level tables ----------
rank_levels <- c("Phylum","Class","Order","Family","Genus","Species")
rank_low <- tolower(rank_levels)

for (i in seq_along(rank_levels)) {
  lv <- rank_levels[i]
  lv_low <- rank_low[i]

  ab <- agg_by(lv, df %>% select(all_of(sample_cols)))
  write_tsv_xls(ab, file.path(OUT_DIR, paste0("allsample.", lv_low, ".abundance.xls")))

  rea <- ab
  rea[sample_cols] <- to_relab(as.matrix(ab[sample_cols]))
  write_tsv_xls(rea, file.path(OUT_DIR, paste0("allsample.", lv_low, ".reabundance.xls")))
}

message("[DONE] allsample tables written to: ", OUT_DIR)
