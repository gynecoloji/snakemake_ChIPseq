#!/usr/bin/env Rscript
# Peak annotation + GO enrichment for one ChIP-seq IP sample.
#
# ChIPseeker annotates each peak to its nearest gene / genomic feature (using the
# UCSC hg38 knownGene TxDb) and clusterProfiler runs GO (BP) enrichment on the
# target genes. Driven by Snakemake's `script:` directive (reads the `snakemake`
# object). Empty / peak-less inputs (e.g. an IgG control) degrade gracefully.

suppressMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(ggplot2)
  library(GenomicRanges)
})

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peaks_file <- snakemake@input[["peaks"]]
anno_out   <- snakemake@output[["anno"]]
feat_out   <- snakemake@output[["featdist"]]
tss_out    <- snakemake@output[["disttss"]]
go_out     <- snakemake@output[["go"]]
name       <- snakemake@params[["name"]]

dir.create(dirname(anno_out), recursive = TRUE, showWarnings = FALSE)

blank <- function(path) ggsave(path, ggplot() + theme_void(), width = 6, height = 4)
empty_go <- function(path)
  write.table(data.frame(ID = character(), Description = character(), p.adjust = numeric()),
              path, sep = "\t", quote = FALSE, row.names = FALSE)

peaks <- tryCatch(readPeakFile(peaks_file), error = function(e) GRanges())

if (length(peaks) == 0) {
  message("No peaks for ", name, " — writing empty annotation/GO outputs.")
  write.table(data.frame(), anno_out, sep = "\t", quote = FALSE, row.names = FALSE)
  empty_go(go_out); blank(feat_out); blank(tss_out)
  quit(save = "no")
}

anno <- annotatePeak(peaks, TxDb = txdb, tssRegion = c(-3000, 3000),
                     annoDb = "org.Hs.eg.db", level = "gene", verbose = FALSE)
df <- as.data.frame(anno)
write.table(df, anno_out, sep = "\t", quote = FALSE, row.names = FALSE)

png(feat_out, width = 1400, height = 700, res = 150)
print(plotAnnoBar(anno) + ggtitle(paste0(name, " — genomic feature distribution")))
dev.off()

png(tss_out, width = 1400, height = 700, res = 150)
print(plotDistToTSS(anno) + ggtitle(paste0(name, " — distance to TSS")))
dev.off()

genes <- unique(stats::na.omit(df$geneId))
ego <- tryCatch(
  enrichGO(gene = genes, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "ENTREZID",
           pAdjustMethod = "BH", qvalueCutoff = 0.2, readable = TRUE),
  error = function(e) NULL)

if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
  write.table(as.data.frame(ego), go_out, sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  message("No significant GO terms for ", name, ".")
  empty_go(go_out)
}
