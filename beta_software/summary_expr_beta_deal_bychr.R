#!/usr/bin/R
args <- commandArgs()
chr <- as.numeric(args[6])
i <- chr
colnames <- read.table("/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/GTEX/v6/weights/Whole_Blood/expr_beta_colnames.txt",header=F,stringsAsFactors=F,sep="\t")
file=paste0("/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/GTEX/v6/weights/Whole_Blood/expr_beta_chr",chr,".txt")
tmp <- read.table(file,header=T,stringsAsFactors=F,sep="\t")
gene_adding <- colnames(tmp)[!colnames(tmp) %in% c("rsID","ref","alt")]

result_df <- data.frame(rsID = tmp$rsID, ref=tmp$ref, alt=tmp$alt)
dat <- matrix(0,nrow=length(tmp$rsID),ncol=(length(colnames$V1)-3))

dat[,match(gene_adding,colnames[-c(1:3)])] <- tmp[,!colnames(tmp) %in% c("rsID","ref","alt")] 
result_df <- cbind(result_df,data.frame(dat))

write.table(result_df,file=paste0("/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/GTEX/v6/weights/Whole_Blood/expr_beta_chr",chr,"_withZero.txt"),col.names=T,row.names=F,sep="\t",quote=F)
