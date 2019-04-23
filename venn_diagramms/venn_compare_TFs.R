setwd("/Users/marasteiger/Documents/Uni/Functional Genomics/peak_calling/ESC/")
cMyc_ESC_path <- "/Users/marasteiger/Documents/Uni/Functional Genomics/peak_calling/ESC/ESC_cMyc_ChIP_seq_peaks.narrowPeak"
cMyc_ESC <- toGRanges(cMyc_ESC_path, format="BED", header=FALSE)

setwd("/Users/marasteiger/Documents/Uni/Functional Genomics/peak_calling/48h/ChipSeq")
cMyc_48h_path <- "/Users/marasteiger/Documents/Uni/Functional Genomics/peak_calling/48h/ChipSeq/48h_cMyc_OSKM_ChIP_seq_peaks.narrowPeak"
cMyc_48h <- toGRanges(cMyc_48h_path, format="BED", header=FALSE)

cMyc_ESC.import <- import(cMyc_ESC_path, format="BED")
cMyc_48h.import <- import(cMyc_48h_path, format="BED")

ol <- findOverlapsOfPeaks(cMyc_ESC,cMyc_48h,connectedPeaks = "merge")
pdf(file="/Users/marasteiger/Documents/Uni/Functional Genomics/peak_calling/venn_TFs_cMyc.pdf")
makeVennDiagram(ol,main="cMyc: ESC vs. 48hMEF",fill = c( "green", "blue"))
dev.off()