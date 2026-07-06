# =============================================================================
# 00_metadata_MR_GR_fullstats_P_SMD.R  (REWRITTEN v4)
# Baseline metadata summary (GR vs MR) with P values, SMD, and Test_method.
#
# v4 update (per user):
# - For Motor phenotype (2 TD; 1 PIGD; 0 indeterminate), display levels as:
#     indeterminate / PIGD / TD
#   (ordered as 0,1,2 respectively), instead of 0/1/2.
#
# User specifications:
# - Exclude "Name" and the grouping column "Levodopa responsiveness (1 MR; 0 GR)" from summary.
# - Treat as Categorical:
#     * "Levodopa (1 Yes; 0 No)"
#     * "Motor phenotype (2 TD; 1 PIGD; 0 indeterminate)"
# - Treat as Count:
#     * "Levodopa/benserazide (tablets/d)"
#     * "Levodopa/carbidopa (tablets/d)"
#
# Default tests:
# - Continuous / Count / Ordinal: Wilcoxon rank-sum test (Mann–Whitney U)
# - Categorical: Fisher's exact test (2×2 or 2×K)
#
# Outputs (in OUT_DIR):
# - Table1_full_stats_P_SMD.csv              (full long table)
# - Table1_main_continuous_ordinal_count.csv (one line per variable)
# - Table1_main_categorical.csv              (level-wise)
# - Imbalance_rank_by_SMD.csv
# - Variable_type_audit.csv
# Stable copy:
# - metadata_full_stats_P_SMD.csv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

## =========================
## 0) Parameters
## =========================
META_FILE <- "metadata.csv"

# Candidate columns that indicate the MR/GR grouping
MR_CANDIDATES <- c(
  "Levodopa responsiveness (1 MR; 0 GR)",
  "MR", "mr",
  "Group", "group",
  "Responsiveness", "response"
)

# Columns to exclude from summary statistics (exact match)
EXCLUDE_EXACT <- c(
  "Name", "name",
  "Sample", "sample", "SampleID", "sampleid",
  "Subject", "subject",
  "ID", "id"
)

# Force ordinal variables by exact name (edit as needed)
ORDINAL_EXACT <- c("H-Y stage", "H-Y Stage", "HY stage", "HY Stage")

# User-specified key variables
V_LEVODOPA_YN <- "Levodopa (1 Yes; 0 No)"
V_MOTOR_PHENO <- "Motor phenotype (2 TD; 1 PIGD; 0 indeterminate)"

# Force categorical variables by exact name
FORCE_CATEGORICAL_EXACT <- c(V_LEVODOPA_YN, V_MOTOR_PHENO)

# Force count variables by exact name
COUNT_EXACT <- c(
  "Levodopa/benserazide (tablets/d)",
  "Levodopa/carbidopa (tablets/d)"
)

# Motor phenotype label mapping (0/1/2 -> labels)
MOTOR_PHENO_CODE_ORDER <- c("0","1","2")
MOTOR_PHENO_LABELS     <- c("indeterminate","PIGD","TD")
names(MOTOR_PHENO_LABELS) <- MOTOR_PHENO_CODE_ORDER

# Output directory
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR <- paste0("Table1_MR_GR_", ts)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(META_FILE))
dt <- fread(META_FILE)

## =========================
## 1) Helpers
## =========================
pick_col <- function(cands, names_vec) {
  hit <- cands[cands %in% names_vec]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

# Detect a "codebook" variable name like:
# "Motor phenotype (2 TD; 1 PIGD; 0 indeterminate)" or "Sex（1 Male; 0 Female）"
# We use presence of digits + semicolon (ASCII or Chinese) + any bracket.
is_codebook_name <- function(v) {
  has_digit <- grepl("[0-9]", v)
  has_semi  <- grepl("[;；]", v)
  has_brkt  <- grepl("[\\(\\（\\[]", v) && grepl("[\\)\\）\\]]", v)
  has_digit && has_semi && has_brkt
}

fmt_mean_sd <- function(x) sprintf("%.2f (%.2f)", mean(x), sd(x))
fmt_median_iqr <- function(x) {
  q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
  sprintf("%.2f [%.2f, %.2f]", q[2], q[1], q[3])
}
fmt_min_max <- function(x) sprintf("%.2f-%.2f", min(x), max(x))
fmt_n_pct <- function(n, N) sprintf("%d (%.1f%%)", n, 100*n/N)

# SMD for continuous-like variables
smd_cont <- function(x, g) {
  x0 <- x[g==0]; x1 <- x[g==1]
  x0 <- x0[!is.na(x0)]; x1 <- x1[!is.na(x1)]
  if (length(x0) < 2 || length(x1) < 2) return(NA_real_)
  m0 <- mean(x0); m1 <- mean(x1)
  s0 <- sd(x0);   s1 <- sd(x1)
  sp <- sqrt((s0^2 + s1^2)/2)
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (m1 - m0)/sp
}

# SMD for categorical variables: compute binary SMD for each level and take max(|SMD|)
smd_cat_max <- function(x, g) {
  x <- as.character(x)
  lv <- sort(unique(x[!is.na(x)]))
  if (length(lv) <= 1) return(NA_real_)
  out <- c()
  for (k in lv) {
    p1 <- mean(x[g==1] == k, na.rm = TRUE)
    p0 <- mean(x[g==0] == k, na.rm = TRUE)
    pbar <- (p1 + p0)/2
    denom <- sqrt(pbar*(1-pbar))
    if (!is.finite(denom) || denom == 0) next
    out <- c(out, (p1 - p0)/denom)
  }
  if (length(out)==0) return(NA_real_)
  max(abs(out))
}

# P values
p_wilcox <- function(x, g) {
  x0 <- x[g==0]; x1 <- x[g==1]
  suppressWarnings(wilcox.test(x1, x0)$p.value)
}
p_fisher_table <- function(tab) suppressWarnings(fisher.test(tab)$p.value)

is_binary_numeric <- function(x) {
  ux <- sort(unique(x[!is.na(x)]))
  length(ux)==2 && all(ux %in% c(0,1))
}

## =========================
## 2) Determine grouping variable (MR vs GR)
## =========================
mr_col <- pick_col(MR_CANDIDATES, names(dt))
if (is.na(mr_col)) stop("Cannot find MR/GR grouping column. Please check metadata.csv column names.")

if (is.numeric(dt[[mr_col]]) || is.integer(dt[[mr_col]])) {
  dt[, MR := as.integer(get(mr_col))]
} else {
  tmp <- tolower(as.character(dt[[mr_col]]))
  dt[, MR := fifelse(tmp %in% c("mr", "moderate", "moderate responder", "mod"), 1L,
                     fifelse(tmp %in% c("gr", "good", "good responder"), 0L, NA_integer_))]
}

dt <- dt[MR %in% c(0,1)]
N0 <- sum(dt$MR==0); N1 <- sum(dt$MR==1)

## =========================
## 3) Variable list (exclude ID/Name + group columns)
## =========================
exclude_cols <- unique(c(
  mr_col, "MR",
  EXCLUDE_EXACT,
  intersect(MR_CANDIDATES, names(dt))
))
vars <- setdiff(names(dt), exclude_cols)

## =========================
## 4) Summarize each variable
## =========================
rows <- list()
audit <- list()

for (v in vars) {
  x <- dt[[v]]
  g <- dt$MR

  miss0 <- sum(is.na(x[g==0]))
  miss1 <- sum(is.na(x[g==1]))

  is_ordinal <- v %in% ORDINAL_EXACT
  is_force_cat <- v %in% FORCE_CATEGORICAL_EXACT
  is_count <- v %in% COUNT_EXACT
  codebook <- is_codebook_name(v)

  inferred_type <- NA_character_
  if (is_force_cat) {
    inferred_type <- "Categorical"
  } else if (is_count) {
    inferred_type <- "Count"
  } else if (is_ordinal) {
    inferred_type <- "Ordinal"
  } else if (is.numeric(x) && !is_binary_numeric(x) && !codebook) {
    inferred_type <- "Continuous"
  } else {
    inferred_type <- "Categorical"
  }

  # Audit record (test method filled below)
  ux <- unique(x[!is.na(x)])
  audit_rec <- data.table(
    Variable = v,
    InferredType = inferred_type,
    is_numeric = is.numeric(x),
    is_binary_numeric = is.numeric(x) && is_binary_numeric(x),
    is_codebook_name = codebook,
    n_unique = length(ux)
  )

  if (inferred_type %in% c("Continuous","Count")) {
    x0 <- x[g==0]; x1 <- x[g==1]
    x0 <- x0[!is.na(x0)]; x1 <- x1[!is.na(x1)]
    if (length(x0)<2 || length(x1)<2) next

    test_method <- "Wilcoxon rank-sum test"

    r <- data.table(
      Variable = v,
      Type = inferred_type,
      Test_method = test_method,
      GR_n = N0, MR_n = N1,
      GR_missing = miss0, MR_missing = miss1,
      GR_mean_SD = fmt_mean_sd(x0),
      MR_mean_SD = fmt_mean_sd(x1),
      GR_median_IQR = fmt_median_iqr(x0),
      MR_median_IQR = fmt_median_iqr(x1),
      GR_min_max = fmt_min_max(x0),
      MR_min_max = fmt_min_max(x1),
      P_value = p_wilcox(x, g),
      SMD = smd_cont(x, g)
    )
    rows[[length(rows)+1]] <- r

    audit_rec[, Test_method := test_method]
    audit[[length(audit)+1]] <- audit_rec

  } else if (inferred_type == "Ordinal") {
    x0 <- x[g==0]; x1 <- x[g==1]
    x0 <- x0[!is.na(x0)]; x1 <- x1[!is.na(x1)]
    if (length(x0)<2 || length(x1)<2) next

    test_method <- "Wilcoxon rank-sum test"

    r <- data.table(
      Variable = v,
      Type = "Ordinal",
      Test_method = test_method,
      GR_n = N0, MR_n = N1,
      GR_missing = miss0, MR_missing = miss1,
      GR_mean_SD = NA_character_,
      MR_mean_SD = NA_character_,
      GR_median_IQR = fmt_median_iqr(x0),
      MR_median_IQR = fmt_median_iqr(x1),
      GR_min_max = fmt_min_max(x0),
      MR_min_max = fmt_min_max(x1),
      P_value = p_wilcox(x, g),
      SMD = smd_cont(x, g)
    )
    rows[[length(rows)+1]] <- r

    audit_rec[, Test_method := test_method]
    audit[[length(audit)+1]] <- audit_rec

  } else {
    # =========================
    # Categorical
    # =========================
    xx <- as.character(x)

    # For output: enforce motor phenotype labels and order
    if (v == V_MOTOR_PHENO) {
      lv_code_out <- MOTOR_PHENO_CODE_ORDER
      lv_label_out <- unname(MOTOR_PHENO_LABELS[lv_code_out])

      # For test: use observed codes only to avoid all-zero columns
      lv_code_test <- sort(unique(xx[!is.na(xx)]))
      lv_code_test <- lv_code_test[lv_code_test %in% lv_code_out]

      if (length(lv_code_test) < 2) {
        pval <- NA_real_
        K <- length(lv_code_test)
      } else {
        tab_test <- table(g, factor(xx, levels = lv_code_test), useNA="no")
        pval <- p_fisher_table(tab_test)
        K <- length(lv_code_test)
      }

      test_method <- if (length(lv_code_test) >= 2) sprintf("Fisher's exact test (2×%d)", K) else "Fisher's exact test"
      smd <- smd_cat_max(xx, g)

      for (i in seq_along(lv_code_out)) {
        k_code <- lv_code_out[i]
        k_label <- lv_label_out[i]

        n0 <- sum(xx[g==0] == k_code, na.rm = TRUE)
        n1 <- sum(xx[g==1] == k_code, na.rm = TRUE)

        r <- data.table(
          Variable = v,
          Level = k_label,
          Type = "Categorical",
          Test_method = test_method,
          GR_n = N0, MR_n = N1,
          GR_missing = miss0, MR_missing = miss1,
          GR_n_pct = fmt_n_pct(n0, N0),
          MR_n_pct = fmt_n_pct(n1, N1),
          P_value = pval,
          SMD = smd
        )
        rows[[length(rows)+1]] <- r
      }

      audit_rec[, Test_method := test_method]
      audit[[length(audit)+1]] <- audit_rec

    } else {
      # Generic categorical variable
      lv_out <- sort(unique(xx[!is.na(xx)]))
      if (length(lv_out) == 0) next

      lv_test <- lv_out
      if (length(lv_test) < 2) {
        pval <- NA_real_
        K <- length(lv_test)
      } else {
        tab_test <- table(g, factor(xx, levels = lv_test), useNA="no")
        pval <- p_fisher_table(tab_test)
        K <- length(lv_test)
      }

      test_method <- if (length(lv_test) >= 2) sprintf("Fisher's exact test (2×%d)", K) else "Fisher's exact test"
      smd  <- smd_cat_max(xx, g)

      for (k in lv_out) {
        n0 <- sum(xx[g==0] == k, na.rm = TRUE)
        n1 <- sum(xx[g==1] == k, na.rm = TRUE)

        r <- data.table(
          Variable = v,
          Level = k,
          Type = "Categorical",
          Test_method = test_method,
          GR_n = N0, MR_n = N1,
          GR_missing = miss0, MR_missing = miss1,
          GR_n_pct = fmt_n_pct(n0, N0),
          MR_n_pct = fmt_n_pct(n1, N1),
          P_value = pval,
          SMD = smd
        )
        rows[[length(rows)+1]] <- r
      }

      audit_rec[, Test_method := test_method]
      audit[[length(audit)+1]] <- audit_rec
    }
  }
}

tab1 <- rbindlist(rows, fill=TRUE)

# Format P and SMD
tab1[, P_value_fmt := ifelse(is.na(P_value), NA_character_,
                             ifelse(P_value < 0.001, "<0.001", sprintf("%.3f", P_value)))]
tab1[, SMD_fmt := ifelse(is.na(SMD), NA_character_, sprintf("%.3f", SMD))]

setorder(tab1, Variable, Type)

audit_dt <- rbindlist(audit, fill=TRUE)
setorder(audit_dt, Variable)
fwrite(audit_dt, file.path(OUT_DIR, "Variable_type_audit.csv"))

## =========================
## 5) Output: full table + split tables
## =========================
fwrite(tab1, file.path(OUT_DIR, "Table1_full_stats_P_SMD.csv"))

# Main table: one line per variable for continuous/ordinal/count
main_num <- tab1[Type %in% c("Continuous","Count","Ordinal"),
                 .(Variable, Type, Test_method,
                   GR = fifelse(Type %in% c("Continuous","Count"), GR_mean_SD, GR_median_IQR),
                   MR = fifelse(Type %in% c("Continuous","Count"), MR_mean_SD, MR_median_IQR),
                   P_value = P_value_fmt,
                   SMD = SMD_fmt)]

main_cat <- tab1[Type=="Categorical",
                 .(Variable, Level, Type, Test_method,
                   GR = GR_n_pct,
                   MR = MR_n_pct,
                   P_value = P_value_fmt,
                   SMD = SMD_fmt)]

fwrite(main_num, file.path(OUT_DIR, "Table1_main_continuous_ordinal_count.csv"))
fwrite(main_cat, file.path(OUT_DIR, "Table1_main_categorical.csv"))

# Imbalance ranking by |SMD|
imb <- unique(tab1[, .(Variable, Type, Test_method, SMD=as.numeric(SMD_fmt), P_value=P_value_fmt)])
imb <- imb[!is.na(SMD)]
imb[, absSMD := abs(SMD)]
setorder(imb, -absSMD)
fwrite(imb, file.path(OUT_DIR, "Imbalance_rank_by_SMD.csv"))

# Stable copy for downstream scripts / tables
file.copy(file.path(OUT_DIR, "Table1_full_stats_P_SMD.csv"),
          "metadata_full_stats_P_SMD.csv",
          overwrite = TRUE)

message("✅ Table 1 recomputed and saved to: ", OUT_DIR)
message("Stable copy written to: metadata_full_stats_P_SMD.csv")
message("Motor phenotype levels displayed as: indeterminate / PIGD / TD.")
