

library(bamsignals)
library(stringr)


### library biomaRt for the genome annotation, we have to choose the corresponding version of the mm9 genome ##
library(biomaRt)
ensembl67 <- useMart(host='may2012.archive.ensembl.org',
                     biomart='ENSEMBL_MART_ENSEMBL',
                     dataset='mmusculus_gene_ensembl')
### now get all (protein coding) genes
prot.gene <- getBM(attributes=c("ensembl_gene_id", "gene_biotype", "strand",
                                "chromosome_name", "start_position", "end_position")
                   ,mart = ensembl67, filters = "biotype", values = "protein_coding")
### now the data frame prot.gene has to be converted into GenomicRanges object. The chromosome names and the strand has to be changed before, as an example, two simple functions....
convert.strand <- function(strand.col){
  xx <-strand.col
  xx <- as.factor(xx)
  levels(xx)[levels(xx)=="-1"] <- "-"
  levels(xx)[levels(xx)=="1"] <- "+"
  levels(xx)[!levels(xx)%in%c("-","+")] <- "*"
  return(xx)
}
###
correct.chr <- function(genes){
  genes <- genes[genes$chromosome_name%in%c(1:19,"X","Y"),]
  genes$chromosome_name <- paste0("chr",genes$chromosome_name)
  return(genes)
}

library(GenomicRanges)

prot.gene$strand <- convert.strand(prot.gene$strand)
prot.gene <- correct.chr(prot.gene)
prot.gene.gr <- makeGRangesFromDataFrame(df=prot.gene, start.field="start_position", end.field="end_position", keep.extra.columns=TRUE,starts.in.df.are.0based=FALSE) #ensembl:1-based

## get the promoters with promoters function
upstr = 2000
downstr = 500
prot.gene.prom <- promoters(prot.gene.gr, upstream = upstr, downstream = downstr)

#counts <- bamCount("ESC_H3K4me3_ChIP-Seq.bam",prot.gene.prom, verbose = FALSE )

files <- c("/project/functional-genomics/2019/group1/mapping/48hMEF/Filtered_ChipSeq/48h_H3.3_ChIP-Seq.bam",
           "/project/functional-genomics/2019/group1/mapping/48hMEF/Filtered_ChipSeq/48h_H3K27ac_ChIP-Seq.bam",
           "/project/functional-genomics/2019/group1/mapping/48hMEF/Filtered_ChipSeq/48h_H3K27me3_ChIP-Seq.bam",
           "/project/functional-genomics/2019/group1/mapping/48hMEF/Filtered_ChipSeq/48h_H3K36me3_ChIP-Seq.bam",
           "/project/functional-genomics/2019/group1/mapping/48hMEF/Filtered_ChipSeq/48h_H3K4me1_ChIP-Seq.bam",
           "/project/functional-genomics/2019/group1/mapping/48hMEF/Filtered_ChipSeq/48h_H3K4me2_ChIP-Seq.bam",
           "/project/functional-genomics/2019/group1/mapping/48hMEF/Filtered_ChipSeq/48h_H3K4me3_ChIP-Seq.bam",
           "/project/functional-genomics/2019/group1/mapping/48hMEF/Filtered_ChipSeq/48h_H3K79me2_ChIP-Seq.bam",
           "/project/functional-genomics/2019/group1/mapping/48hMEF/Filtered_ChipSeq/48h_H3K9ac_ChIP-Seq.bam",
           "/project/functional-genomics/2019/group1/mapping/48hMEF/Filtered_ChipSeq/48h_H3K9me3_ChIP-Seq.bam",
           "/project/functional-genomics/2019/group1/mapping/48hMEF/Filtered_ChipSeq/48h_H3_ChIP-Seq.bam",
           "/project/functional-genomics/2019/group1/mapping/48hMEF/Filtered_ChipSeq/48h_Hdac1_ChIP-Seq.bam",
           "/project/functional-genomics/2019/group1/mapping/48hMEF/Filtered_ChipSeq/48h_Inputnative_MNase_ChIP-Seq.bam")


# files <- c("/project/functional-genomics/2019/group1/mapping/ESC/Filtered_ChipSeq/ESC_H3.3_ChIP-Seq.bam",
#            "/project/functional-genomics/2019/group1/mapping/ESC/Filtered_ChipSeq/ESC_H3_ChIP-Seq.bam",
#            "/project/functional-genomics/2019/group1/mapping/ESC/Filtered_ChipSeq/ESC_H3K4me1_ChIP-Seq.bam",
#            "/project/functional-genomics/2019/group1/mapping/ESC/Filtered_ChipSeq/ESC_H3K4me2_ChIP-Seq.bam",
#            "/project/functional-genomics/2019/group1/mapping/ESC/Filtered_ChipSeq/ESC_H3K4me3_ChIP-Seq.bam",
#            "/project/functional-genomics/2019/group1/mapping/ESC/Filtered_ChipSeq/ESC_H3K9ac_ChIP-Seq.bam",
#            "/project/functional-genomics/2019/group1/mapping/ESC/Filtered_ChipSeq/ESC_H3K9me3_ChIP-Seq.bam",
#            "/project/functional-genomics/2019/group1/mapping/ESC/Filtered_ChipSeq/ESC_H3K27ac_ChIP-Seq.bam",
#            "/project/functional-genomics/2019/group1/mapping/ESC/Filtered_ChipSeq/ESC_H3K27me3_ChIP-Seq.bam",
#            "/project/functional-genomics/2019/group1/mapping/ESC/Filtered_ChipSeq/ESC_H3K36me3_ChIP-Seq.bam",
#            "/project/functional-genomics/2019/group1/mapping/ESC/Filtered_ChipSeq/ESC_H3K79me2_ChIP-Seq.bam",
#            "/project/functional-genomics/2019/group1/mapping/ESC/Filtered_ChipSeq/ESC_Hdac1_ChIP-Seq.bam",
#            "/project/functional-genomics/2019/group1/mapping/ESC/Filtered_ChipSeq/ESC_Inputnative_MNase_ChIP-Seq.bam")

setwd("/project/functional-genomics/2019/group1/linear_regression")

# read counts from files by using the bamCounts function and promoter regions from mart 
counts <- c()
for (i in 1:length(files)){
  counts <- cbind(counts, bamCount(files[i] , prot.gene.prom, verbose = FALSE ))
}

# add histone modification names to columns and write raw counts to table
names <- (str_extract(files,"H.*_"))
names <- substr(names,1,nchar(names)-1)
names[length(names)] <- "MNase"
colnames(counts) <- names
write.table(counts,"counts_48h.txt",row.names=F, quote=F)

# extract control counts to seperate vector
cntrl <- ncol(counts)
control <- counts[,cntrl]
counts <- counts[,-cntrl]

# calculate median of counts per promoter region
med <- c()
for (p in 1:ncol(counts)){
  med <- append(med,median((counts[,p]+1)/(control+1)))
}

for (h in 1:ncol(counts)){
  counts[,h] <- ((counts[,h]+1)/(control+1))*(1/med[h]) 
}

# normalize counts 
#for (p in 1:nrow(counts)){
#  for (h in 1:ncol(counts)){
#    counts[p,h] <- ((counts[p,h]+1)/(control[p]+1)) * (1/med[p])
#  }
#}
counts <- log2(counts)


norm <- counts

#scale counts in copy called norm and write to new file with normalized counts to be used by linear regression script 
for (p in 1:ncol(norm)){
  # d <- mean(norm[p,])
  # norm[p,] <- scale(norm[p,], scale=F, center=d-1)
  # norm[p,] <- scale(norm[p,], scale=T, center=F)
  norm[,p] <- scale(norm[,p], scale = T, center = F)
  norm[,p] <- scale(norm[,p], scale = F, center = mean(norm[,p])-1)
}

write.table(norm,"norm_48h_new2.txt",row.names=F, quote=F, col.names=names[1:ncol(norm)])


