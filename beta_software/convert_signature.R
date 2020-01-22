#!/usr/bin/R

signature <- read.table("/ysm-gpfs/pi/zhao-data/mc2792/sc_gwas/deconv/cybersort/LM22.txt",header=T,stringsAsFactors=F,sep="\t")
cell_types <- colnames(signature)[-1]
genes <- signature[,1]
dat <- t(signature[,-1])
result <- cbind(data.frame(cell.type=cell_types),dat)
colnames(result) <- c("cell.type",genes)
write.table(result,file="/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/signature/LM22.txt",col.names=T,row.names=F,sep="\t",quote=F)
