#!/usr/bin/R
library("genio")
library(glmnet)

k=10
alpha = 0.5
args <- commandArgs()
chr <- args[6]
gene <- args[7]
tissue <- args[8]
gene_i <- gsub("\\.[0-9]+","",gene)
cis_snp_dir <- paste0("/ysm-gpfs/project/wl382/GTEx_v8/genotype/cis_loc/",chr,"/",gene,"/",gene)
expr_dir <- paste0("/gpfs/loomis/scratch60/fas/radev/zy92/GTEX/adjusted_expr1/",chr,"/",gene_i,"/",tissue,".adj_expr")

geno <- read_plink(cis_snp_dir)
expr <- read.table(expr_dir,header=F,stringsAsFactors=F)

bed <- geno$X # snp by sample
bim <- data.frame(geno$bim)
fam <- data.frame(geno$fam)

# deal with the sample id in expr to match them with that in geno
a <- strsplit(expr[,1],"-")
sample_expr <- sapply(a,function(X) paste0(X[1],"-",X[2]))
inter_fam <- intersect(sample_expr,fam[,1])


# replace na with rowMeans
row_means <-  apply(bed,1,function(X) X[is.na(X)] <- mean(na.omit(X)))
inds <- which(is.na(bed),arr.ind=T)
bed[inds] <- row_means[inds[,1]]

#get the data for intersected samples
expr_1 <- expr[match(inter_fam,sample_expr),]
bed_1 <- bed[,match(inter_fam,fam[,1])]
fam_1 <- fam[match(inter_fam,fam[,1]),]

fit <- cv.glmnet(t(bed_1),expr_1[,2],nfolds=k,alpha=alpha,keep=T,parallel=F) 
fit.df <- data.frame(fit$cvm,fit$lambda,1:length(fit$cvm)) 
best.lam <- fit.df[which.min(fit.df[,1]),]
cvm.best = best.lam[,1]
lambda.best = best.lam[,2]
nrow.best = best.lam[,3] ##position of best lambda in cv.glmnet output
      
ret <- as.data.frame(fit$glmnet.fit$beta[,nrow.best]) # get betas from best lambda
ret[ret == 0.0] <- NA
bestbetas = as.vector(ret[which(!is.na(ret)),]) # vector of non-zero betas
names(bestbetas) = bim[!is.na(ret),2]
if(length(bestbetas)==0){
  next()
}

# for bim file (A1:non-effective; A2: effective)
result <- bim[!is.na(ret),]
result$weight <- bestbetas
write.table(result,file=paste0("/ysm-gpfs/scratch60/wl382/GTEx/v8/models/",chr,"/",gene,"/",tissue,"_weight.txt"),col.names=T,row.names=F,sep="\t",quote=F)

pred.mat <- fit$fit.preval[,nrow.best]
res <- summary(lm(expr_1[,2]~pred.mat))
rsq <- res$r.squared
pval <- res$coef[2,4]
result_df <- data.frame(gene=gene,tissue=tissue,rsq=rsq,pval=pval)
write.table(result_df,file=paste0("/ysm-gpfs/scratch60/wl382/GTEx/v8/models/",chr,"/",gene,"/",tissue,"_rsq.txt"),col.names=T,row.names=F,sep="\t",quote=F)
