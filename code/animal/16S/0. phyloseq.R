## 00_import_phyloseq.R
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(phyloseq)
})

FEATURE_TSV <- "featuretable_no_singletons.txt"
OUT_RDS     <- "phyloseq_16S.rds"

stopifnot(file.exists(FEATURE_TSV))

# 1) 读取：跳过第一行 "# Constructed from biom file"
df <- readr::read_tsv(FEATURE_TSV, skip = 1, show_col_types = FALSE)

# 2) 统一列名
colnames(df)[1] <- "FeatureID"
stopifnot("taxonomy" %in% names(df))

sample_cols <- setdiff(names(df), c("FeatureID", "taxonomy"))

# 3) 构建 sample metadata（由样本名首字母映射分组）
group_of <- function(x){
  if (startsWith(x, "A")) return("GR")
  if (startsWith(x, "B")) return("MR")
  if (startsWith(x, "C")) return("MR+Phage cocktail")
  if (startsWith(x, "D")) return("MR+AFMT")
  return(NA_character_)
}

meta <- tibble(
  sample = sample_cols,
  group  = vapply(sample_cols, group_of, character(1))
) %>% mutate(group = factor(group, levels = c("GR","MR","MR+Phage cocktail","MR+AFMT")))

if (any(is.na(meta$group))) stop("Some samples cannot be mapped to group. Check sample names.")

# 4) OTU 矩阵
otu_mat <- df %>%
  select(FeatureID, all_of(sample_cols)) %>%
  column_to_rownames("FeatureID") %>%
  as.matrix()

# 5) taxonomy 解析：k__/p__/c__/o__/f__/g__/s__
# 注意你的 taxonomy 形式是 "k__Bacteria; p__Firmicutes; ..."
tax_raw <- df$taxonomy %>% as.character()

split_tax <- str_split(tax_raw, pattern = ";\\s*", simplify = TRUE)
# 保守取前 7 列（k p c o f g s）；不足会填空
if (ncol(split_tax) < 7) {
  split_tax <- cbind(split_tax, matrix("", nrow(split_tax), 7 - ncol(split_tax)))
} else {
  split_tax <- split_tax[, 1:7, drop = FALSE]
}

colnames(split_tax) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

strip_prefix <- function(x){
  x <- str_trim(x)
  x <- str_replace(x, "^[a-z]__", "")
  x <- na_if(x, "")
  # 把 unclassified / uncultured 等统一为 NA（也可以保留字符串，看你偏好）
  ifelse(is.na(x) | str_detect(tolower(x), "unclassified|uncultured"), NA, x)
}

tax_mat <- apply(split_tax, 2, strip_prefix)
rownames(tax_mat) <- df$FeatureID

# 6) phyloseq 对象
ps <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(meta %>% column_to_rownames("sample"))
)

saveRDS(ps, OUT_RDS)
cat("[OK] Saved:", OUT_RDS, "\n")
