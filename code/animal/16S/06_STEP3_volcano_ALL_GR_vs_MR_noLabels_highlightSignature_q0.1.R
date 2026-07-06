## STEP3_volcano_yQ_GR_vs_MR_highlightSignature_q0.1_NO_RING.R
## Volcano (ALL taxa) for 6 levels: GR vs MR
## Y = -log10(q)  (back to previous version)
## No taxa labels; highlight fingerprint taxa by larger SAME-COLOR solid points (NO black ring)

suppressPackageStartupMessages({
  library(ggplot2)
})

# ---------------- SETTINGS ----------------
SIG_DIR <- "16S_Signature_q0.1_6Levels"     # 这里放 ALL_* 和 SIGNATURE_* 的csv
OUTDIR  <- "16S_Volcano_ALL_GR_vs_MR_yQ_q0.1_6Levels_NO_RING"
dir.create(OUTDIR, showWarnings = FALSE)

LEVELS <- c("phylum","class","order","family","genus","species")

Q_CUT <- 0.10
LOG2FC_CUT <- 1.0          # “log2FC默认”：常用阈值±1；不想要竖线可改成 0 或注释掉 vline

# colors
COL_DOWN <- "#2C7BB6"
COL_UP   <- "#D7191C"
COL_NS   <- "grey40"

BASE_SIZE <- 2.2
SIG_SIZE  <- 3.4           # 指纹点高亮大小（实心点，无描边）
ALPHA_ALL <- 0.85
ALPHA_SIG <- 0.95

LOGFILE <- file.path(OUTDIR, "volcano_log.txt")
cat("", file = LOGFILE)  # reset log
log_msg <- function(...) cat(paste0(..., "\n"), file = LOGFILE, append = TRUE)
# -----------------------------------------

sanitize_names <- function(nm) {
  nm[is.na(nm)] <- ""
  nm <- trimws(nm)
  blank <- which(nm == "")
  if (length(blank) > 0) nm[blank] <- paste0("V", blank)
  make.names(nm, unique = TRUE)
}

safe_read_csv <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- sanitize_names(names(df))
  df
}

get_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) NA_character_ else hit[1]
}

find_all_file <- function(level) {
  f1 <- file.path(SIG_DIR, paste0("ALL_", level, "_GR_vs_MR_CLR_Wilcox_BH.csv"))
  if (file.exists(f1)) return(f1)
  cand <- list.files(SIG_DIR, pattern = paste0("^ALL_", level, ".*\\.csv$"), full.names = TRUE)
  if (length(cand) > 0) return(cand[1])
  NA_character_
}

read_signature_taxa <- function(level) {
  f <- file.path(SIG_DIR, paste0("SIGNATURE_q0.1_", level, "_GR_vs_MR.csv"))
  if (!file.exists(f)) return(character(0))
  sig <- tryCatch(safe_read_csv(f), error = function(e) NULL)
  if (is.null(sig)) return(character(0))
  if (!("taxon" %in% names(sig))) return(character(0))
  unique(as.character(sig$taxon))
}

make_one <- function(level) {
  all_file <- find_all_file(level)
  if (is.na(all_file)) {
    log_msg("[SKIP] ", level, ": cannot find ALL file under ", SIG_DIR)
    return(invisible(NULL))
  }

  df <- safe_read_csv(all_file)

  # columns (Step1 standard / tolerant)
  col_taxon <- get_col(df, c("taxon"))
  col_fc    <- get_col(df, c("log2FC_rel","log2fc_rel","log2FC","log2fc"))
  col_q     <- get_col(df, c("q_value","q","fdr","BH_q","qvalue"))

  if (any(is.na(c(col_taxon, col_fc, col_q)))) {
    log_msg("[FAIL] ", level, ": required columns missing in ", basename(all_file))
    log_msg("       available cols: ", paste(names(df), collapse = ", "))
    stop("Missing required columns at level=", level)
  }

  taxon <- as.character(df[[col_taxon]])
  log2fc <- suppressWarnings(as.numeric(df[[col_fc]]))
  q <- suppressWarnings(as.numeric(df[[col_q]]))

  ok <- !is.na(taxon) & taxon != "" & !is.na(log2fc) & !is.na(q)
  taxon <- taxon[ok]; log2fc <- log2fc[ok]; q <- q[ok]

  q_plot <- pmax(q, 1e-300)
  neglog10q <- -log10(q_plot)

  pass_q <- q < Q_CUT
  direction <- rep("Not significant", length(q))
  direction[pass_q & log2fc >=  LOG2FC_CUT] <- "Up (q<0.1)"
  direction[pass_q & log2fc <= -LOG2FC_CUT] <- "Down (q<0.1)"
  direction <- factor(direction, levels = c("Down (q<0.1)","Not significant","Up (q<0.1)"))

  sig_taxa <- read_signature_taxa(level)
  is_signature <- taxon %in% sig_taxa

  dat <- data.frame(
    taxon = taxon,
    log2FC = log2fc,
    neglog10q = neglog10q,
    q_value = q,
    direction = direction,
    is_signature = is_signature,
    stringsAsFactors = FALSE
  )

  p <- ggplot(dat, aes(x = log2FC, y = neglog10q)) +
    geom_hline(yintercept = -log10(Q_CUT), linetype = "dashed", linewidth = 0.5, color = "grey40") +
    geom_vline(xintercept = c(-LOG2FC_CUT, LOG2FC_CUT), linetype = "dashed", linewidth = 0.5, color = "grey40") +

    # all points
    geom_point(aes(color = direction), size = BASE_SIZE, alpha = ALPHA_ALL, shape = 16) +

    # highlight signature points: SAME color, bigger solid points (NO ring, NO stroke)
    geom_point(
      data = dat[dat$is_signature, , drop = FALSE],
      aes(color = direction),
      size = SIG_SIZE, alpha = ALPHA_SIG, shape = 16
    ) +

    scale_color_manual(values = c(
      "Down (q<0.1)" = COL_DOWN,
      "Not significant" = COL_NS,
      "Up (q<0.1)" = COL_UP
    )) +
    theme_classic(base_size = 13) +
    theme(legend.title = element_blank()) +
    labs(x = "log2FC (GR vs MR)", y = expression(-log[10](q)), title = NULL)

  out_png <- file.path(OUTDIR, paste0("Volcano_ALL_", level, "_GR_vs_MR_yQ_q0.1_noLabels_NO_RING.png"))
  out_pdf <- file.path(OUTDIR, paste0("Volcano_ALL_", level, "_GR_vs_MR_yQ_q0.1_noLabels_NO_RING.pdf"))

  ggsave(out_png, p, width = 7.6, height = 5.6, dpi = 300)
  ggsave(out_pdf, p, width = 7.6, height = 5.6)

  log_msg("[OK] ", level,
          " n=", nrow(dat),
          " sig(q&|fc|)=", sum(dat$q_value < Q_CUT & abs(dat$log2FC) >= LOG2FC_CUT),
          " signature=", sum(dat$is_signature),
          " -> ", basename(out_png))
}

options(error = NULL)  # 防止递归错误处理导致“严重bug”
for (lv in LEVELS) {
  tryCatch(
    make_one(lv),
    error = function(e) log_msg("[ERROR] ", lv, " :: ", conditionMessage(e))
  )
}

message("[DONE] Output folder: ", OUTDIR, " | log: ", LOGFILE)
