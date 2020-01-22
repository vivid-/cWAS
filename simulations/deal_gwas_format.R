#!/usr/bin/R

args <- commandArgs()
gwasf <- args[6]
bimf <- args[7]
outputf <- args[8]

gwas <- read.table(gwasf,header=T,stringsAsFactors=F)
bim <- read.table(bimf,header=F,stringsAsFactors=F)

rsID <- gwas[,"SNP"]
Z <- gwas[,"STAT"]
N <- gwas[,"NMISS"]
alt <- gwas[,"A1"]
ref <- bim[match(rsID,bim[,2]),6]
result_df <- data.frame(rsID=rsID,ref=ref,alt=alt,Z=Z,N=N)
write.table(result_df,file=outputf,col.names=T,row.names=F,sep="\t",quote=F)
