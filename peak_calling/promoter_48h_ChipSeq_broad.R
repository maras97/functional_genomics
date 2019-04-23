library(biomaRt)
library(GenomicFeatures)

setwd("/Users/marasteiger/Documents/Uni/Functional Genomics/data")
mart = useMart('ENSEMBL_MART_ENSEMBL',dataset='mmusculus_gene_ensembl', host="may2012.archive.ensembl.org")

txdb = makeTxDbFromGFF("Mus_musculus.NCBIM37.67.gtf", format="gtf")
PR <- promoters(txdb, upstream=2000)

setwd("/Users/marasteiger/Documents/Uni/Functional Genomics/promoter_peaks/")
tabelle <- read.table("48h_ChIP.txt",header=F,stringsAsFactors = F)
names <- tabelle[,1]
names<- unlist(lapply(names,paste,"_peaks.narrowPeak"))
names <- sub(" ","",names)
counts <- c()
length <- c()

broad <- c("48h_H3K9me3_ChIP-Seq_peaks.narrowPeak",
           "48h_H3K27me3_ChIP-Seq_peaks.narrowPeak",
           "48h_H3K36me3_ChIP-Seq_peaks.narrowPeak",
           "48h_H3K79me2_ChIP-Seq_peaks.narrowPeak" )
for (n in 1:length(names)){
  if (names[n] %in% broad){
    names[n] <- sub("_peaks.narrowPeak","_peaks.broadPeak", names[n])
  }
}

setwd("/Users/marasteiger/Documents/Uni/Functional Genomics/peak_calling/48h/ChipSeq")

for (filename in names){
  if (!grepl("Input",filename)){
    file_path <- paste("/Users/marasteiger/Documents/Uni/Functional Genomics/peak_calling/48h/ChipSeq/",filename,sep="")
    file <- toGRanges(file_path, format="BED", header=FALSE)
    seqlevels(file) <- sub("chr", "", seqlevels(file))
    seqlevels(file) <- sub("M", "MT", seqlevels(file))
    ol <- findOverlaps(file,PR,ignore.strand=F,type="within")
    peaks_promoters <- length(unique(queryHits(ol)))
    length <- append(length,queryLength(ol))
    counts <- append(counts,peaks_promoters)
    print(counts)
  }
  else if (grepl("Input",filename)){
    counts <- append(counts,0)
    length <- append(length,0)
  }
}
tabelle <- cbind(tabelle,length,counts)




bartable = rbind(tabelle[,4],tabelle[,5],nrow=2)
cols <- sub("ESC_","",names)
cols <- sub("48h_","",cols)
cols <- sub("_ChIP-Seq_peaks.narrowPeak","",cols)
cols <- sub("_ChIP_seq_peaks.narrowPeak","",cols)
cols <- sub("_ChIP-seq_peaks.narrowPeak","",cols)
cols <- sub("_ChIP-Seq_peaks.broadPeak","",cols)
cols <- sub("Inputnative_","",cols)
colnames(bartable) <- cols
par(mar = c(6, 4, 2, 2))
barplot(bartable, beside = TRUE,las=2,cex.names=0.65)
legend("topright",  fill=c("gray27","grey"), c("Peaks", "Anteil in Promotoren")  )

setwd("/Users/marasteiger/Documents/Uni/Functional Genomics/promoter_peaks/")
write.table(tabelle,"48h_ChIP_broad.txt")
