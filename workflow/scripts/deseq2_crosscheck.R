#!/usr/bin/env Rscript
# Independent cross-check of the pydeseq2 differential-expression result using
# Bioconductor DESeq2. Reads the pseudobulk counts + sample metadata exported by
# notebooks/04_differential_expression.ipynb, refits the paired ~ donor + group model,
# and reports the overlap of significant genes with the Python result.
#
# Usage:
#   Rscript deseq2_crosscheck.R pseudobulk_counts.csv pseudobulk_coldata.csv out_overlap.csv
#
# Install (once):
#   BiocManager::install("DESeq2")

suppressMessages({
  library(DESeq2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: deseq2_crosscheck.R <counts.csv> <coldata.csv> <out.csv>")
}
counts_path  <- args[1]
coldata_path <- args[2]
out_path     <- args[3]

FDR    <- 0.05
LOG2FC <- 1.0

# Load pseudobulk (samples x genes) and sample metadata
counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
coldata <- read.csv(coldata_path, row.names = 1, check.names = FALSE)
mat <- round(t(as.matrix(counts)))
mode(mat) <- "integer"
coldata <- coldata[colnames(mat), , drop = FALSE]
stopifnot(all(c("donor", "group") %in% colnames(coldata)))
stopifnot(identical(colnames(mat), rownames(coldata)))
coldata$donor <- factor(coldata$donor)

# Detect the contrast from the group levels so the same script serves the dorsolateral
# primary and the anteroposterior robustness pseudobulk. The reference level is the negative
# arm (ventromedial or anterior), matching pydeseq2's contrast order.
levels_present <- sort(unique(as.character(coldata$group)))
ref_level <- if ("ventromedial" %in% levels_present) "ventromedial" else "anterior"
stopifnot(ref_level %in% levels_present)
pos_level <- setdiff(levels_present, ref_level)[1]
py_file <- if (ref_level == "ventromedial") {
  "de_dorsolateral_vs_ventromedial_deseq2.csv"
} else {
  "de_posterior_vs_anterior_deseq2.csv"
}
py_path <- file.path(dirname(out_path), py_file)
coldata$group <- relevel(factor(coldata$group), ref = ref_level)

# Paired design, matching the Python primary model
dds <- DESeqDataSetFromMatrix(countData = mat, colData = coldata,
                              design = ~ donor + group)
dds <- DESeq(dds, quiet = TRUE)
res <- as.data.frame(results(dds, contrast = c("group", pos_level, ref_level)))
res$gene <- rownames(res)
sig_r <- subset(res, !is.na(padj) & padj < FDR & abs(log2FoldChange) > LOG2FC)$gene

# Compare with pydeseq2
overlap_n <- NA_integer_
jaccard   <- NA_real_
sig_py <- character(0)
if (file.exists(py_path)) {
  py <- read.csv(py_path)
  gcol <- intersect(c("gene_identifier", "gene", "X", "names", "index"), colnames(py))[1]
  pcol <- intersect(c("padj", "pvalue_adj", "qvalue"), colnames(py))[1]
  fcol <- intersect(c("log2FoldChange", "log2fc", "lfc"), colnames(py))[1]
  if (!is.na(gcol) && !is.na(pcol) && !is.na(fcol)) {
    sig_py <- py[!is.na(py[[pcol]]) & py[[pcol]] < FDR &
                   abs(py[[fcol]]) > LOG2FC, gcol]
    inter <- length(intersect(sig_r, sig_py))
    uni   <- length(union(sig_r, sig_py))
    overlap_n <- inter
    jaccard   <- if (uni > 0) inter / uni else NA_real_
  }
} else {
  message("Python cross-check result not found at ", py_path, "; reporting R-only counts.")
}

summary_df <- data.frame(
  metric = c("n_sig_deseq2_R", "n_sig_pydeseq2", "n_overlap", "jaccard"),
  value  = c(length(sig_r), length(sig_py), overlap_n, round(jaccard, 4))
)
write.csv(summary_df, out_path, row.names = FALSE)
cat(sprintf("DESeq2(R) sig=%d  pydeseq2 sig=%d  overlap=%s  jaccard=%s\n",
            length(sig_r), length(sig_py),
            as.character(overlap_n), as.character(round(jaccard, 4))))
