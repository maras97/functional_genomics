setwd("/Users/marasteiger/Documents/Uni/Functional Genomics/peak_calling/48h/ATACSeq")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ChIPpeakAnno", version = "3.8")
require("ChIPpeakAnno")

MEF <- "/Users/marasteiger/Documents/Uni/Functional Genomics/peak_calling/48h/ChipSeq/48h_H3K27ac_ChIP-Seq_summits.bed"
MEF <- toGRanges(MEF, format="BED", header=FALSE)

ESC <- "/Users/marasteiger/Documents/Uni/Functional Genomics/peak_calling/ESC/ChipSeq/ESC_H3K27ac_ChIP-Seq_summits.bed"
ESC <- toGRanges(ESC, format="BED", header=FALSE)

library(rtracklayer)
MEF.import <- import(MEF, format="BED")
ESC.import <- import(ESC, format="BED")
identical(start(MEF), start(MEF.import))
identical(start(ESC), start(ESC.import))
#gff <- "/Users/marasteiger/Documents/Uni/Functional Genomics/peak_calling/48h/ATACSeq/48h_KpMX_ATAC-seq_summits.bed"
#gr2 <- toGRanges(gff, format="GFF", header=FALSE, skip=3)
ol <- findOverlapsOfPeaks(MEF, ESC)
pdf(file="/Users/marasteiger/Documents/Uni/Functional Genomics/peak_calling/ESC/ChipSeq/venn_H3K27ac_ChIP-Seq.pdf")
makeVennDiagram(ol,main="H3K27ac_ChIP-Seq")
dev.off()

