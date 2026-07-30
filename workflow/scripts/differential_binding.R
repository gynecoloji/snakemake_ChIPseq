#!/usr/bin/env Rscript
# Differential binding for one contrast (condition A vs B) via DESeq2 on the
# consensus fragment-count matrix (results/consensus/consensus_counts.txt).
#
# log2FoldChange > 0 => higher binding in condition A. Significant peaks are
# annotated to genes with ChIPseeker. Driven by Snakemake's `script:` directive.
# NOTE: contrasts are pre-filtered in common.smk to those with >=2 replicates per
# condition (DESeq2 needs replicates for dispersion), so this only runs on
# adequately-replicated designs.

suppressMessages({
  library(DESeq2)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(GenomicRanges)
})

counts_file <- snakemake@input[["counts"]]
peaks_bed   <- snakemake@input[["peaks"]]
res_out     <- snakemake@output[["results"]]
sig_out     <- snakemake@output[["sig"]]
ma_out      <- snakemake@output[["ma"]]
volcano_out <- snakemake@output[["volcano"]]
pca_out     <- snakemake@output[["pca"]]
cond_a <- snakemake@params[["cond_a"]]
cond_b <- snakemake@params[["cond_b"]]
bams_a <- unlist(snakemake@params[["bams_a"]])
bams_b <- unlist(snakemake@params[["bams_b"]])
name   <- snakemake@params[["name"]]

dir.create(dirname(res_out), recursive = TRUE, showWarnings = FALSE)

# featureCounts matrix: 6 annotation columns then one column per BAM (the column
# name is the BAM path passed to featureCounts, e.g. results/blacklist_filtered/<s>.nobl.bam).
fc <- read.delim(counts_file, comment.char = "#", check.names = FALSE)
rownames(fc) <- fc$Geneid
ann_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
mat <- as.matrix(fc[, !(colnames(fc) %in% ann_cols), drop = FALSE])

sel <- c(bams_a, bams_b)
sel <- sel[sel %in% colnames(mat)]
if (length(sel) < 2) stop("differential_binding: <2 matching count columns for ", name)
m <- round(mat[, sel, drop = FALSE])
cond <- factor(ifelse(sel %in% bams_a, cond_a, cond_b), levels = c(cond_b, cond_a))
coldata <- data.frame(row.names = sel, condition = cond)

dds <- DESeqDataSetFromMatrix(countData = m, colData = coldata, design = ~condition)
dds <- tryCatch(DESeq(dds), error = function(e) {
  message("DESeq() failed (", conditionMessage(e), "); retrying with local fit + gene-est dispersion")
  DESeq(dds, fitType = "local")
})
res <- results(dds, contrast = c("condition", cond_a, cond_b))
resdf <- as.data.frame(res)
resdf$peak <- rownames(resdf)
resdf <- resdf[, c("peak", setdiff(colnames(resdf), "peak"))]
write.table(resdf, res_out, sep = "\t", quote = FALSE, row.names = FALSE)

sig <- subset(resdf, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
if (nrow(sig) > 0) {
  peaks <- tryCatch(readPeakFile(peaks_bed), error = function(e) GRanges())
  if (length(peaks) > 0 && !is.null(peaks$V4)) {
    anno <- as.data.frame(annotatePeak(
      peaks[peaks$V4 %in% sig$peak], TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
      tssRegion = c(-3000, 3000), annoDb = "org.Hs.eg.db", level = "gene", verbose = FALSE))
    key <- if ("V4" %in% colnames(anno)) "V4" else NULL
    if (!is.null(key)) sig <- merge(sig, anno[, c(key, "annotation", "SYMBOL", "distanceToTSS")],
                                    by.x = "peak", by.y = key, all.x = TRUE)
  }
}
write.table(sig, sig_out, sep = "\t", quote = FALSE, row.names = FALSE)

png(ma_out, width = 1200, height = 900, res = 150)
plotMA(res, main = paste0(name, "  (", cond_a, " vs ", cond_b, ")"), ylim = c(-5, 5))
dev.off()

resdf$sig <- with(resdf, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
png(volcano_out, width = 1100, height = 950, res = 150)
print(ggplot(resdf, aes(log2FoldChange, -log10(pmax(padj, 1e-300)), color = sig)) +
        geom_point(alpha = 0.5, size = 0.8) +
        scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "#c62828"), guide = "none") +
        labs(x = "log2 fold change", y = "-log10 adjusted p", title = name) +
        theme_bw())
dev.off()

vsd <- tryCatch(varianceStabilizingTransformation(dds, blind = TRUE),
                error = function(e) normTransform(dds))
png(pca_out, width = 1000, height = 900, res = 150)
print(plotPCA(vsd, intgroup = "condition") + ggtitle(name) + theme_bw())
dev.off()
