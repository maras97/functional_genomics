setwd("/Users/marasteiger/Documents/Uni/Functional Genomics/peak_calling/48h/ChipSeq")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ChIPpeakAnno", version = "3.8")
require("ChIPpeakAnno")

oct4_path <- "/Users/marasteiger/Documents/Uni/Functional Genomics/peak_calling/48h/48h_Oct4_OSKM_ChIP_seq_peaks.narrowPeak"
oct4 <- toGRanges(oct4_path, format="BED", header=FALSE)
Sox2_path <- "/Users/marasteiger/Documents/Uni/Functional Genomics/peak_calling/48h/ChipSeq/48h_Sox2_OSKM_ChIP_seq_peaks.narrowPeak"
Sox2 <- toGRanges(Sox2_path, format="BED", header=FALSE)
Klf4_path <- "/Users/marasteiger/Documents/Uni/Functional Genomics/peak_calling/48h/ChipSeq/48h_Klf4_OSKM_ChIP_seq_peaks.narrowPeak"
Klf4 <- toGRanges(Klf4_path, format="BED", header=FALSE)
cMyc_path <- "/Users/marasteiger/Documents/Uni/Functional Genomics/peak_calling/48h/ChipSeq/48h_cMyc_OSKM_ChIP_seq_peaks.narrowPeak"
cMyc <- toGRanges(cMyc_path, format="BED", header=FALSE)



library(rtracklayer)
oct4.import <- import(oct4_path, format="BED")
Sox2.import <- import(Sox2_path, format="BED")
Klf4.import <- import(Klf4_path, format="BED")
cMyc.import <- import(cMyc_path, format="BED")

#gff <- "/Users/marasteiger/Documents/Uni/Functional Genomics/peak_calling/48h/ATACSeq/48h_KpMX_ATAC-seq_summits.bed"
#gr2 <- toGRanges(gff, format="GFF", header=FALSE, skip=3)
ol <- findOverlapsOfPeaks(oct4, Sox2,Klf4,cMyc)
pdf(file="/Users/marasteiger/Documents/Uni/Functional Genomics/peak_calling/48h/ChipSeq/venn_TFs_48h_ChIP-Seq.pdf")
makeVennDiagram(ol,main="ESC_TFs_ChipSeq",fill = c("darkblue", "darkgreen", "orange", "darkorchid4"))
dev.off()


#pdf(file="/Users/marasteiger/Documents/Uni/Functional Genomics/peak_calling/ESC/pie_TFs_ESC_ChIP-Seq.pdf")
#pie1(table(ol$overlappingPeaks[["oct4///Sox2"]]$overlapFeature))
#dev.off()

overlappingPeaks = ol$OverlappingPeaks

