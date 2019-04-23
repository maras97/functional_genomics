setwd("/project/functional-genomics/2019/group1/mapping/48hMEF/RNASeq")
names <- read.table(file="/project/functional-genomics/2019/group1/rename.txt",fill = T, stringsAsFactors = F, sep="", comment.char = "", check.names = FALSE, header = F)

files <- list.files(pattern="*SRR*", full.names=FALSE, recursive=F)

new_name<-c()
for(i in 1:length(files)){
  a <- files[i]
  start <- regexpr("SRR[0-9]*",a)[1]
  length <- attr(regexpr("SRR[0-9]*",a),"match.length")
  file_id <- substr(a,start,start+length-1)
  row <- which(names[,2]==file_id)
  file_new <- gsub(file_id,names[row,1], files[i])
  new_name <- append(new_name,file_new)
}
print(files)
print(new_name)
file.rename(files,new_name)
