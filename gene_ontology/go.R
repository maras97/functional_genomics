require(GSEABase)
require(Biobase)
library(gsean)
library(WGCNA)
biocLite("DESeq2")
library(DESeq2)
library("vsn")
library(biomaRt)
library(GenomicFeatures)
require(org.Mm.eg.db)

setwd("/Users/marasteiger/Documents/Uni/Functional Genomics/data")
mart = useMart('ENSEMBL_MART_ENSEMBL',dataset='mmusculus_gene_ensembl', host="may2012.archive.ensembl.org")

setwd("/Users/marasteiger/Documents/Uni/Functional Genomics/mapping/")
data <- read.table('fusionPE.txt', header = T, fill=T, stringsAsFactors = F)
data <- data[,c(-5:-9)]
colnames(data)<- sub("X","",colnames(data))
countdata <- as.matrix(data)
condition <- factor(c(rep("48h",4), rep("ESC",2)))
type <- factor(c(rep("paired-end",1),c(rep("single-read",4), rep("paired-end",1))))
coldata <- data.frame(condition,type)
dds <- DESeqDataSetFromMatrix(countdata, colData = coldata, design = ~ condition + type)
dds <- DESeq(dds)
r <- results(dds, alpha = 0.01,pAdjustMethod="BH",contrast = c("condition","ESC","48h"))


setwd("/Users/marasteiger/Documents/Uni/Functional Genomics/gene_ontology/")

keep <- (r$log2FoldChange) < 5 & r$padj < 0.01
sum(keep,na.rm=TRUE)
keep[is.na(keep)] <- FALSE

dds2 <- dds[keep,]

dds2 <- DESeq(dds2)
resultsNames(dds2)
## [1] "conditiontreated"   "conditionuntreated"
res <- results(dds2, contrast = c("condition", "ESC","48h") )
statistic <- res$stat
names(statistic) <- rownames(res)
exprs <- counts(dds2, normalized = TRUE)

# convert gene id
rows <- c()
for (k in 1:length(names(statistic))){
  name <- (getBM(attributes=c('external_gene_id'),filters='ensembl_gene_id',
                 values = names(statistic)[k],mart=mart)$external_gene_id)
  rows <- append(rows,name)
}
names(statistic) <- rows
rownames(exprs) <- rows


require(org.Mm.eg.db)

#gene.id <- AnnotationDbi::select(org.Mm.eg.db, names(statistic),"ENSEMBL", "ENTREZID","GENENAME")

library(gskb)

#load(system.file("extdata","Mus_musculus.gmt", package = "gsean"))

data(mm_GO)

#require(GSA)
#gmt <- GSA.read.gmt("Mus_musculus.gmt")


# GSEA

set.seed(1)
# result.GSEA <- gsean(gmt, statistic, exprs)
# result.GSEA <- gsean(GO_dme, statistic, exprs)

result.GSEA <- gsean(mm_GO, statistic, exprs)

invisible(capture.output(p <- GSEA.barplot(result.GSEA, category = 'pathway',
                                           score = 'NES', top = 50, pvalue = 'padj',
                                           sort = 'padj', numChar = 110) + 
                           theme(plot.margin = margin(10, 10, 10, 50))))
plotly::ggplotly(p)

# GSEA.barplot(result.GSEA, category = 'pathway',
#              score = 'NES', top = 50, pvalue = 'padj',
#              sort = 'padj', numChar = 110) 





