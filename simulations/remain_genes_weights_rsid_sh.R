#!/usr/bin/R

dat <- read.table("remain_WB_genes_weights.txt",header=T,stringsAsFactors=F,sep="\t")
sink("remain_WB_genes_weight_reformat.sh")
for(i in 1:nrow(dat)){
#  chr <- paste0("chr",dat[i,1])
  chr <- dat[i,1]
  gene <- dat[i,2]
  sent <- paste0("Rscript ~/weightRsid_byGene.R Whole_Blood ",chr," ",gene)
  cat(sent,"\n")
}
sink()
