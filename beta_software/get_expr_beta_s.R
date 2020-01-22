#!/usr/bin/R

args <- commandArgs(trailingOnly = TRUE)
# pre calculate the \hat\beta %*% S in each chromosome
out_f <- strsplit(grep('--out*', args, value = TRUE), split = '=')[[1]][[2]]
expr_beta_f <- strsplit(grep('--expr-beta*', args, value = TRUE), split = '=')[[1]][[2]]
S_f <- strsplit(grep('--signature-matrix*', args, value = TRUE), split = '=')[[1]][[2]]


# read expr_beta_f # only for matched genes in the signature matrix
# format: SNP.ID REF ALT GENE1 GENE2 GENE3...
expr_beta <- read.table(expr_beta_f,header=T,stringsAsFactors=F,sep="\t")
# read signatrue matrix in the corresponding tissue
#format: CELL.type GENE1 GENE2 GENE3 ...
signature <- read.table(S_f,header=T,stringsAsFactors=F,sep="\t")

# match gene in expr_beta and signature
inter_gene <- intersect(colnames(expr_beta)[-c(1:3)],colnames(signature)[-1])
# match intersected genes
expr_beta_1 <- expr_beta[,c(1:3,match(inter_gene,colnames(expr_beta)))]
signature_1 <- signature[,c(1,match(inter_gene,colnames(signature)))]
#\hat\beta %*% S
S <- t(signature_1[,-1])
beta <- expr_beta[,-c(1:3)]
result <- data.frame(as.matrix(beta) %*% as.matrix(S))
result$rsID <- expr_beta[,1]
result$ref <- expr_beta[,2]
result$alt <- expr_beta[,3]
colnames(result) <- c(signature[,1],"rsID","ref","alt")
write.table(result,file = out_f,col.names=T,row.names=F,sep="\t",quote=F)
