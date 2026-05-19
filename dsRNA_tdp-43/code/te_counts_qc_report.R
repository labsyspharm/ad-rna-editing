#!/usr/bin/env Rscript
# ============================================================================
# te_counts_qc_report.R
# QC report generator for featureCounts (Rsubread) TE/repeat quantification
# output saved as a .qs file.
#
# Usage:
#   Rscript te_counts_qc_report.R \
#       --fc        /path/to/amp-ad_repeat_counts_raw-fraction_AND_secondary_noMETA.qs \
#       --rmsk      /path/to/repeatmasker_raw.csv \
#       --metadata  /path/to/sample_metadata.csv   (optional)
#       --out       /path/to/qc_report.pdf
#
# Designed for fc_results produced with:
#   countMultiMappingReads=TRUE, fraction=TRUE, primaryOnly=FALSE,
#   allowMultiOverlap=TRUE, isPairedEnd=TRUE
# ============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(qs)
  library(tidyverse)
  library(matrixStats)
  library(pheatmap)
  library(RColorBrewer)
  library(patchwork)
  library(DESeq2)
})

# ---------------------------------------------------------------- CLI parsing
option_list <- list(
  make_option(c("--fc"),       type = "character", help = "Path to fc_results .qs file"),
  make_option(c("--rmsk"),     type = "character", default = NULL,
              help = "Path to RepeatMasker CSV (for class_family aggregation)"),
  make_option(c("--metadata"), type = "character", default = NULL,
              help = "Optional sample metadata CSV (must contain a sample_id column)"),
  make_option(c("--out"),      type = "character", default = "qc_report.pdf",
              help = "Output PDF path [default %default]"),
  make_option(c("--top_n"),    type = "integer",   default = 20,
              help = "Top-N TE families to plot [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
stopifnot("--fc is required" = !is.null(opt$fc))

# -------------------------------------------------------------------- helpers
log_msg <- function(...) message("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)

clean_sample_names <- function(x) {
  x <- basename(x)
  x <- sub("\\.bam$", "", x)
  x <- sub("\\.final$", "", x)
  x
}

safe_pct <- function(num, den) ifelse(den > 0, num / den * 100, NA_real_)

# ---------------------------------------------------------------------- load
log_msg("Loading fc_results from ", opt$fc)
fc <- qread(opt$fc)
stopifnot(all(c("counts", "annotation", "stat") %in% names(fc)))

counts <- fc$counts
colnames(counts) <- clean_sample_names(colnames(counts))

stat <- fc$stat
colnames(stat)[-1] <- clean_sample_names(colnames(stat)[-1])

annot <- as_tibble(fc$annotation)

log_msg("Counts matrix: ", nrow(counts), " features x ", ncol(counts), " samples")

# Optional RepeatMasker for class/family aggregation
rmsk <- NULL
if (!is.null(opt$rmsk) && file.exists(opt$rmsk)) {
  log_msg("Loading RepeatMasker annotation from ", opt$rmsk)
  rmsk <- read_csv(opt$rmsk, show_col_types = FALSE) %>%
    mutate(GeneID = paste(rep_name, row_number(), sep = "_")) %>%
    select(GeneID, rep_name, class_family)
}

# Optional metadata
meta <- NULL
if (!is.null(opt$metadata) && file.exists(opt$metadata)) {
  log_msg("Loading metadata from ", opt$metadata)
  meta <- read_csv(opt$metadata, show_col_types = FALSE)
}

# =============================================================== compute QC
# Per-sample summary
log_msg("Computing per-sample summaries")
lib_size       <- colSums(counts)
detected       <- colSums(counts > 0)
total_reads    <- colSums(stat[,-1])
assigned       <- as.numeric(stat[stat$Status == "Assigned", -1])
unmapped_no_ft <- as.numeric(stat[stat$Status == "Unassigned_NoFeatures", -1])
unmapped_mm    <- as.numeric(stat[stat$Status == "Unassigned_MultiMapping", -1])
unmapped_amb   <- as.numeric(stat[stat$Status == "Unassigned_Ambiguity", -1])

sample_qc <- tibble(
  sample_id       = colnames(counts),
  total_reads     = total_reads,
  assigned        = assigned,
  assignment_pct  = safe_pct(assigned, total_reads),
  lib_size_counts = lib_size,
  detected_loci   = detected,
  no_feature_pct  = safe_pct(unmapped_no_ft, total_reads),
  multimap_pct    = safe_pct(unmapped_mm,    total_reads),
  ambig_pct       = safe_pct(unmapped_amb,   total_reads)
)

# Flag outliers (>3 MAD from median on key metrics)
flag <- function(x) abs(x - median(x, na.rm = TRUE)) > 3 * mad(x, na.rm = TRUE)
sample_qc <- sample_qc %>%
  mutate(flag_libsize  = flag(lib_size_counts),
         flag_assign   = flag(assignment_pct),
         flag_detected = flag(detected_loci),
         is_outlier    = flag_libsize | flag_assign | flag_detected)

n_outliers <- sum(sample_qc$is_outlier, na.rm = TRUE)
log_msg("Flagged ", n_outliers, " potential outlier samples")

# =================================================================== plots
log_msg("Building plots")
theme_set(theme_minimal(base_size = 10) +
            theme(panel.grid.minor = element_blank()))

# 1. Stat: stacked bar of assignment categories
stat_long <- stat %>%
  pivot_longer(-Status, names_to = "sample_id", values_to = "n") %>%
  filter(n > 0) %>%
  group_by(sample_id) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p_stat <- ggplot(stat_long,
                 aes(x = reorder(sample_id, -n, sum), y = pct, fill = Status)) +
  geom_col(width = 1) +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Read-assignment status per sample",
       subtitle = "Each column is one sample. Look for samples whose Assigned fraction differs from the cohort.",
       x = NULL, y = "% of reads") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom")

# 2. Library size distribution
p_lib <- ggplot(sample_qc, aes(lib_size_counts)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  geom_vline(xintercept = median(sample_qc$lib_size_counts), linetype = 2, color = "red") +
  labs(title = "Per-sample total assigned counts",
       subtitle = "Bimodality usually signals a batch effect.",
       x = "Sum of counts (this is fractional, so it's an effective fragment count)",
       y = "Samples")

# 3. Assignment rate vs library size
p_scatter <- ggplot(sample_qc, aes(total_reads, assignment_pct, color = is_outlier)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "grey30"), name = "Outlier") +
  labs(title = "Assignment rate vs. total reads",
       subtitle = "Samples with low % assigned but normal read count = something is wrong upstream.",
       x = "Total reads in BAM (sum of stat)", y = "% Assigned to repeats")

# 4. Detection rate vs library size
p_detect <- ggplot(sample_qc, aes(lib_size_counts, detected_loci, color = is_outlier)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "grey30"), name = "Outlier") +
  labs(title = "Detected loci vs. library size",
       subtitle = "High library + low detection = signal concentrated on a few loci (bad).",
       x = "Sum of counts", y = "Number of features with count > 0")

# 5. TE class composition (requires rmsk)
p_class <- NULL
p_top   <- NULL
class_long <- NULL
if (!is.null(rmsk)) {
  log_msg("Aggregating to TE class / family")
  rmsk_match <- rmsk %>% filter(GeneID %in% rownames(counts))
  ord <- match(rownames(counts), rmsk_match$GeneID)
  class_vec  <- rmsk_match$class_family[ord]
  family_vec <- rmsk_match$rep_name[ord]

  class_counts <- rowsum(counts, group = class_vec, na.rm = TRUE)
  class_long <- as_tibble(class_counts, rownames = "class_family") %>%
    pivot_longer(-class_family, names_to = "sample_id", values_to = "n") %>%
    group_by(sample_id) %>% mutate(pct = n / sum(n) * 100) %>% ungroup()

  # Keep only top classes for readability
  top_classes <- class_long %>% group_by(class_family) %>%
    summarise(m = mean(pct)) %>% slice_max(m, n = 10) %>% pull(class_family)
  class_long_top <- class_long %>%
    mutate(class_family = ifelse(class_family %in% top_classes, class_family, "Other"))

  p_class <- ggplot(class_long_top,
                    aes(reorder(sample_id, -n, sum), pct, fill = class_family)) +
    geom_col(width = 1) +
    scale_fill_brewer(palette = "Paired") +
    labs(title = "TE class composition per sample",
         subtitle = "Should be stable across samples. Outliers in composition are red flags.",
         x = NULL, y = "% of assigned counts") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "bottom")

  # Top families
  fam_counts <- rowsum(counts, group = family_vec, na.rm = TRUE)
  fam_mean <- sort(rowMeans(fam_counts), decreasing = TRUE)
  top_fams <- names(head(fam_mean, opt$top_n))
  fam_long <- as_tibble(fam_counts[top_fams, ], rownames = "rep_name") %>%
    pivot_longer(-rep_name, names_to = "sample_id", values_to = "n")

  p_top <- ggplot(fam_long, aes(reorder(rep_name, n, median), n)) +
    geom_boxplot(outlier.size = 0.5, fill = "grey80") +
    scale_y_log10() + coord_flip() +
    labs(title = paste0("Top ", opt$top_n, " TE families by mean count"),
         subtitle = "Box = across-sample distribution. Wide spread = candidate for batch confounding.",
         x = NULL, y = "Counts (log10)")
}

# 6. Sample-sample correlation heatmap (subsample features for speed)
log_msg("Computing sample correlation matrix")
set.seed(1)
keep_features <- which(rowSums(counts) > 0)
n_sub <- min(50000, length(keep_features))
sub_idx <- sample(keep_features, n_sub)
log_counts <- log1p(counts[sub_idx, , drop = FALSE])
cor_mat <- cor(log_counts)

cor_long <- as_tibble(cor_mat, rownames = "a") %>%
  pivot_longer(-a, names_to = "b", values_to = "r")

p_cor_hist <- ggplot(cor_long %>% filter(a != b), aes(r)) +
  geom_histogram(bins = 50, fill = "darkorange", color = "white") +
  labs(title = "Sample-sample correlation distribution",
       subtitle = "Tight high-r cluster is good. A long left tail = candidate outliers.",
       x = "Pearson r (log1p counts)", y = "Sample pairs")

# 7. PCA via DESeq2 vst
log_msg("Running VST + PCA")
counts_int <- round(counts)
storage.mode(counts_int) <- "integer"
keep_genes <- rowSums(counts_int) >= 10
counts_int <- counts_int[keep_genes, , drop = FALSE]

col_df <- data.frame(sample_id = colnames(counts_int), row.names = colnames(counts_int))
if (!is.null(meta)) {
  col_df <- col_df %>% left_join(meta, by = "sample_id")
  rownames(col_df) <- col_df$sample_id
}

dds <- DESeqDataSetFromMatrix(counts_int, colData = col_df, design = ~ 1)
vsd <- tryCatch(vst(dds, blind = TRUE),
                error = function(e) { log_msg("vst() failed, falling back to varianceStabilizingTransformation");
                  varianceStabilizingTransformation(dds, blind = TRUE) })

vst_mat <- assay(vsd)
top_var <- order(rowVars(vst_mat), decreasing = TRUE)[seq_len(min(5000, nrow(vst_mat)))]
pca <- prcomp(t(vst_mat[top_var, ]))
pct_var <- round(100 * (pca$sdev^2) / sum(pca$sdev^2), 1)

pca_df <- as_tibble(pca$x[, 1:5]) %>% mutate(sample_id = colnames(vst_mat))
if (!is.null(meta)) pca_df <- pca_df %>% left_join(meta, by = "sample_id")
pca_df <- pca_df %>% left_join(sample_qc %>% select(sample_id, is_outlier), by = "sample_id")

p_pca12 <- ggplot(pca_df, aes(PC1, PC2, color = is_outlier)) +
  geom_point(alpha = 0.85) +
  scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "grey30"), name = "Flagged") +
  labs(title = "PCA (PC1 vs PC2)",
       subtitle = sprintf("PC1 %.1f%%   PC2 %.1f%%   |   top 5000 variable features after VST",
                          pct_var[1], pct_var[2]),
       x = sprintf("PC1 (%.1f%%)", pct_var[1]),
       y = sprintf("PC2 (%.1f%%)", pct_var[2]))

p_pca23 <- ggplot(pca_df, aes(PC2, PC3, color = is_outlier)) +
  geom_point(alpha = 0.85) +
  scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "grey30"), name = "Flagged") +
  labs(title = "PCA (PC2 vs PC3)",
       x = sprintf("PC2 (%.1f%%)", pct_var[2]),
       y = sprintf("PC3 (%.1f%%)", pct_var[3]))

# Additional PCA panels colored by metadata variables (categorical w/ <=10 levels)
meta_pcas <- list()
if (!is.null(meta)) {
  candidate_cols <- setdiff(colnames(meta), "sample_id")
  for (cc in candidate_cols) {
    v <- pca_df[[cc]]
    if (is.null(v)) next
    if (is.numeric(v)) {
      meta_pcas[[cc]] <- ggplot(pca_df, aes(PC1, PC2, color = .data[[cc]])) +
        geom_point(alpha = 0.85) +
        scale_color_viridis_c() +
        labs(title = paste0("PCA colored by ", cc))
    } else if (length(unique(na.omit(v))) <= 10) {
      meta_pcas[[cc]] <- ggplot(pca_df, aes(PC1, PC2, color = factor(.data[[cc]]))) +
        geom_point(alpha = 0.85) +
        labs(title = paste0("PCA colored by ", cc), color = cc)
    }
  }
}

# ================================================================ write PDF
log_msg("Writing report to ", opt$out)
pdf(opt$out, width = 11, height = 8.5)

# Cover
grid::grid.newpage()
grid::grid.text(
  paste0("featureCounts TE Quantification QC Report\n",
         "Input: ", basename(opt$fc), "\n",
         "Features: ", format(nrow(counts), big.mark = ","), "  |  Samples: ", ncol(counts), "\n",
         "Outliers flagged: ", n_outliers, "\n",
         "Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M")),
  gp = grid::gpar(fontsize = 14)
)

print(p_stat)
print(p_lib + p_scatter + plot_layout(ncol = 2))
print(p_detect + p_cor_hist + plot_layout(ncol = 2))
if (!is.null(p_class)) print(p_class)
if (!is.null(p_top))   print(p_top)

# Correlation heatmap (cap at 100 samples for legibility)
if (ncol(cor_mat) <= 100) {
  pheatmap(cor_mat, show_rownames = FALSE, show_colnames = FALSE,
           color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
           main = "Sample-sample correlation heatmap")
} else {
  # Downsample
  set.seed(2); sub_s <- sample(colnames(cor_mat), 100)
  pheatmap(cor_mat[sub_s, sub_s], show_rownames = FALSE, show_colnames = FALSE,
           color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
           main = "Sample correlation heatmap (100-sample subsample)")
}

print(p_pca12 + p_pca23 + plot_layout(ncol = 2))
for (pp in meta_pcas) print(pp)

# Outlier table
if (n_outliers > 0) {
  outliers <- sample_qc %>% filter(is_outlier)
  gridExtra::grid.table(outliers %>% select(sample_id, lib_size_counts, assignment_pct,
                                            detected_loci, flag_libsize, flag_assign, flag_detected))
}

invisible(dev.off())

# Also write the sample_qc tibble as CSV alongside the PDF
csv_out <- sub("\\.pdf$", "_sample_qc.csv", opt$out)
write_csv(sample_qc, csv_out)
log_msg("Wrote per-sample QC table to ", csv_out)
log_msg("Done.")
