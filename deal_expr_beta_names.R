#!/usr/bin/R

args <- commandArgs()
genes <- read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/gtexGene.txt",header=F,stringsAsFactors=F)
expr_beta <- read.table(args[6],header=T,stringsAsFactors=F,sep="\t")

ensembl <- colnames(expr_beta)[-c(1:3)]
gene_ids <- genes[match(ensembl,genes[,7]),4]
gene_ids[is.na(gene_ids)] <- ensembl[is.na(gene_ids)]
colnames(expr_beta) <- c(colnames(expr_beta)[c(1:3)],gene_ids)
write.table(expr_beta,file=args[7],col.names=T,row.names=F,sep="\t",quote=F)
