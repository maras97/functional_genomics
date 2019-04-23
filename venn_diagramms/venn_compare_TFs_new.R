if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ChIPpeakAnno", version = "3.8")
require("ChIPpeakAnno")
library(rtracklayer)
library(eulerr)

setwd("/Users/marasteiger/Documents/Uni/Functional Genomics/peak_calling/ESC/")
Klf4_ESC_path <- "/Users/marasteiger/Documents/Uni/Functional Genomics/peak_calling/ESC/ESC_Klf4_ChIP_seq_peaks.narrowPeak"

Klf4_ESC <- toGRanges(Klf4_ESC_path, format="BED", header=FALSE)

setwd("/Users/marasteiger/Documents/Uni/Functional Genomics/peak_calling/48h/ChipSeq")
Klf4_48h_path <- "/Users/marasteiger/Documents/Uni/Functional Genomics/peak_calling/48h/ChipSeq/48h_Klf4_OSKM_ChIP_seq_peaks.narrowPeak"
Klf4_48h <- toGRanges(Klf4_48h_path, format="BED", header=FALSE)

Klf4_ESC.import <- import(Klf4_ESC_path, format="BED")
Klf4_48h.import <- import(Klf4_48h_path, format="BED")



ol <- findOverlapsOfPeaks(Klf4_ESC,Klf4_48h)
a <- ol[["all.peaks"]][["Klf4_ESC"]]@elementMetadata@nrows
b <- ol[["all.peaks"]][["Klf4_48h"]]@elementMetadata@nrows
ab <- nrow(ol[["overlappingPeaks"]][["Klf4_ESC///Klf4_48h"]])
fit <- euler(c(A = a-ab,  B= b-ab, "A&B" = ab))

pdf(file="/Users/marasteiger/Documents/Uni/Functional Genomics/peak_calling/venn_TFs_Klf4_new.pdf")
plot(fit,labels = c("Klf4_ESC","Klf4_48h",""),quantities=T, fill=c("lightblue","lightgreen"))
dev.off()


