# the script to fit cell-type interactive imputation model
# !/usr/bin/R
library("genio")
library(glmnet)

k=10
alpha = 0.5
args <- commandArgs()
chr <- args[6]
gene <- args[7]
tissue <- args[8]
tissue_frac <- args[9]
set.seed(1)
assayed_frac_dir = paste0("/gpfs/ysm/pi/zhao-data/wl382/cWAS/data/assayed_fractions/",tissue_frac,"/")
#assayed_frac_f <- paste0(assayed_frac_dir,tissue_frac,"_CiberSort_frac.txt")
assayed_frac_f <- paste0(assayed_frac_dir,"assayed_genes_expr_frac.txt")
assayed_frac <- read.table(assayed_frac_f,header=T,stringsAsFactors=F,sep="\t")
cols_2remove <- c("P.value","Pearson.Correlation","RMSE")
assayed_frac <- assayed_frac[,which(!(colnames(assayed_frac) %in% cols_2remove))]
# deal with the fam
#frac_ids <- assayed_frac[,1]
frac_ids <- assayed_frac$FID
frac_ids_tmp <- strsplit(frac_ids,"_")
frac_ids_tmp_1 <- sapply(frac_ids_tmp,"[[",1)
frac_ids_tmp_2 <- sapply(frac_ids_tmp,"[[",2)
frac_ids_tmp_str <- paste0(frac_ids_tmp_1,"-",frac_ids_tmp_2)
#assayed_frac[,1] <- frac_ids_tmp_str
assayed_frac$FID <- frac_ids_tmp_str

gene_i <- gsub("\\.[0-9]+","",gene)
cis_snp_dir <- paste0("/ysm-gpfs/project/wl382/GTEx_v8/genotype/cis_loc/",chr,"/",gene,"/",gene)

expr_dir <- paste0("/gpfs/loomis/scratch60/zhao/wl382/GTEx_V8/adjusted_expr/",chr,"/")

gene_i <- list.files(expr_dir,pattern=gene_i)
expr_dir <- paste0(expr_dir,gene_i,"/",tissue,".adj_expr")
geno <- read_plink(cis_snp_dir)
expr <- read.table(expr_dir,header=F,stringsAsFactors=F)

bed <- geno$X # snp by sample
bim <- data.frame(geno$bim)
fam <- data.frame(geno$fam)

# deal with the sample id in expr to match them with that in geno
a <- strsplit(expr[,1],"-")
sample_expr <- sapply(a,function(X) paste0(X[1],"-",X[2]))
inter_fam <- intersect(sample_expr,fam[,1])
inter_fam <- intersect(inter_fam,assayed_frac$FID)
#inter_fam <- intersect(inter_fam,assayed_frac[,1])


# replace na with rowMeans
row_means <-  apply(bed,1,function(X) X[is.na(X)] <- mean(na.omit(X)))
inds <- which(is.na(bed),arr.ind=T)
bed[inds] <- row_means[inds[,1]]

#get the data for intersected samples
expr_1 <- expr[match(inter_fam,sample_expr),]
bed_1 <- bed[,match(inter_fam,fam[,1])]
fam_1 <- fam[match(inter_fam,fam[,1]),]

#assayed_frac_1 <- assayed_frac[match(inter_fam,assayed_frac[,1]),]
assayed_frac_1 <- assayed_frac[match(inter_fam,assayed_frac$FID),]

expr_1 <- expr[match(inter_fam,sample_expr),]

geno  <- t(bed_1)
n_cell <- ncol(assayed_frac_1)-1
geno_frac <- data.frame()
for(i in 1:n_cell){
	#tmp <- assayed_frac_1[,i+1] * geno
	tmp <- assayed_frac_1[,i] * geno
	if(i == 1){
		geno_frac <- tmp
	}else{
		geno_frac <- cbind(geno_frac,tmp)
	}
}

X <- geno_frac
#X <- cbind(geno,geno_frac)
X <- as.matrix(X)

fit <- cv.glmnet(X,scale(expr_1[,2],center=T,scale=F),nfolds=k,alpha=alpha,keep=T,parallel=F) 
fit.df <- data.frame(fit$cvm,fit$lambda,1:length(fit$cvm)) 
best.lam <- fit.df[which.min(fit.df[,1]),]
cvm.best = best.lam[,1]
lambda.best = best.lam[,2]
nrow.best = best.lam[,3] ##position of best lambda in cv.glmnet output
#
pred.mat <- fit$fit.preval[,nrow.best]
res <- summary(lm(expr_1[,2]~pred.mat))
rsq <- res$r.squared
pval <- res$coef[2,4]

result_df_model1 <- data.frame(gene=gene,tissue=tissue,rsq=rsq,pval=pval,model="cell_interactive_QTL_nonConstraint_only")

############################################################################

assayed_frac_dir = paste0("/gpfs/ysm/pi/zhao-data/wl382/cWAS/data/assayed_fractions/",tissue_frac,"/")

assayed_frac_f <- paste0(assayed_frac_dir,tissue_frac,"_CiberSort_frac.txt")
assayed_frac <- read.table(assayed_frac_f,header=T,stringsAsFactors=F,sep="\t")
cols_2remove <- c("P.value","Pearson.Correlation","RMSE")
assayed_frac <- assayed_frac[,which(!(colnames(assayed_frac) %in% cols_2remove))]
## deal with the fam
frac_ids <- assayed_frac[,1]
frac_ids_tmp <- strsplit(frac_ids,"_")
frac_ids_tmp_1 <- sapply(frac_ids_tmp,"[[",1)
frac_ids_tmp_2 <- sapply(frac_ids_tmp,"[[",2)
frac_ids_tmp_str <- paste0(frac_ids_tmp_1,"-",frac_ids_tmp_2)
assayed_frac[,1] <- frac_ids_tmp_str

inter_fam <- intersect(sample_expr,fam[,1])
inter_fam <- intersect(inter_fam,assayed_frac[,1])
#
##get the data for intersected samples
expr_1 <- expr[match(inter_fam,sample_expr),]
bed_1 <- bed[,match(inter_fam,fam[,1])]
fam_1 <- fam[match(inter_fam,fam[,1]),]

assayed_frac_1 <- assayed_frac[match(inter_fam,assayed_frac[,1]),]

expr_1 <- expr[match(inter_fam,sample_expr),]

geno  <- t(bed_1)
n_cell <- ncol(assayed_frac_1)-1
geno_frac <- data.frame()
for(i in 1:n_cell){
	tmp <- assayed_frac_1[,i+1] * geno
	if(i == 1){
		geno_frac <- tmp
	}else{
		geno_frac <- cbind(geno_frac,tmp)
	}
}

X <- geno_frac
#X <- cbind(geno,geno_frac)
X <- as.matrix(X)
fit <- cv.glmnet(X,scale(expr_1[,2],center=T,scale=F),nfolds=k,alpha=alpha,keep=T,parallel=F) 
fit.df <- data.frame(fit$cvm,fit$lambda,1:length(fit$cvm)) 
best.lam <- fit.df[which.min(fit.df[,1]),]
cvm.best = best.lam[,1]
lambda.best = best.lam[,2]
nrow.best = best.lam[,3] ##position of best lambda in cv.glmnet output

pred.mat <- fit$fit.preval[,nrow.best]
res <- summary(lm(expr_1[,2]~pred.mat))
rsq <- res$r.squared
pval <- res$coef[2,4]

result_df_model0 <- data.frame(gene=gene,tissue=tissue,rsq=rsq,pval=pval,model="cell_interactive_QTL_CiberSort_only")


############################################################################


# the second model is impute gene expression based on imputed cell type fraction
imputed_frac_dir <- paste0("/gpfs/ysm/pi/zhao-data/wl382/cWAS/data/imputed_fractions/",tissue_frac,"/")
imputed_frac_f <- paste0(imputed_frac_dir,"predicted_genes_expr_frac.txt")
imputed_frac <- read.table(imputed_frac_f,header=T,stringsAsFactors=F,sep="\t")
imputed_frac_1 <- imputed_frac[match(inter_fam,imputed_frac$FID),]

#sig_dir <- "/gpfs/ysm/pi/zhao-data/wl382/cWAS/data/signature_denoised/"
#expr_files <- list.files(sig_dir,pattern=tissue)
#fetal_type = gsub(paste0("HCL_",tissue),"",tissue_frac)
#fetal_type <- gsub("_","",tissue_frac)
#if(fetal_type=="fetal"){
#  fetal_type="Fetal"
#}else{
#  fetal_type="Adult"
#}
#expr_f <- expr_files[grep(fetal_type,expr_files)]
#expr_f <- expr_f[grep("_saverX.txt",expr_f)]
#expr_f <- paste0(sig_dir,expr_f)
#sig_expr <- read.table(expr_f,header=T,stringsAsFactors=F,sep="\t")
#gene_loc <- read.table("/ysm-gpfs/project/wl382/GTEx_v8/gene_names_noVersion_locs_v8.txt",header=T,stringsAsFactors=F,sep="\t")
#col_names <- colnames(sig_expr)
#test <- gene_loc[match(col_names[-1],gene_loc[,3]),2]



#X_1 <- cbind(geno,imputed_frac_1[,which(colnames(imputed_frac_1)!="FID")])
#X_1 <- as.matrix(X_1)
X_1 <- as.matrix(imputed_frac_1[,which(colnames(imputed_frac_1)!="FID")])
X_new <- cbind(geno,X_1)
fit_1 <- cv.glmnet(X_new,scale(expr_1[,2],center=T,scale=F),nfolds=k,alpha=alpha,keep=T,parallel=F) 
fit.df <- data.frame(fit_1$cvm,fit_1$lambda,1:length(fit_1$cvm)) 
best.lam <- fit.df[which.min(fit.df[,1]),]
cvm.best = best.lam[,1]
lambda.best = best.lam[,2]
nrow.best = best.lam[,3] ##position of best lambda in cv.glmnet output

pred.mat <- fit_1$fit.preval[,nrow.best]
res <- summary(lm(expr_1[,2]~pred.mat))
rsq <- res$r.squared
pval <- res$coef[2,4]
result_df_model2 <- data.frame(gene=gene,tissue=tissue,rsq=rsq,pval=pval,model="GRP_withGenetics_only")

result_df_total <- rbind(result_df_model1,result_df_model0)
result_df_total <- rbind(result_df_total,result_df_model2)
#result_df_total <- result_df_model2

result_dir <- paste0("/gpfs/ysm/pi/zhao-data/wl382/cWAS/data/fit_expr/",chr,"/",gene,"/")
dir.create(result_dir,recursive=TRUE,showWarnings =FALSE)
#result_file <- paste0(result_dir,tissue,"_",tissue_frac,"_fitted_expr_r2_noAdditionalGeno.txt")
#write.table(result_df_total,file=result_file,col.names=T,row.names=F,sep="\t",quote=F)


############################################################################
# the third model is to impute gene expression based on imputed cell type fraction and signature gene expression
library(coefplot)
X_new <- cbind(X,X_1)
X_new <- cbind(geno,X_new)
fit_2 <- cv.glmnet(X_new,scale(expr_1[,2],center=T,scale=F),nfolds=k,alpha=alpha,keep=T,parallel=F) 
fit.df <- data.frame(fit_2$cvm,fit_2$lambda,1:length(fit_2$cvm)) 
best.lam <- fit.df[which.min(fit.df[,1]),]
cvm.best = best.lam[,1]
lambda.best = best.lam[,2]
nrow.best = best.lam[,3] 

pred.mat <- fit_1$fit.preval[,nrow.best]
res <- summary(lm(expr_1[,2]~pred.mat))
rsq <- res$r.squared
pval <- res$coef[2,4]
result_df_model3 <- data.frame(gene=gene,tissue=tissue,rsq=rsq,pval=pval,model="GRP_ieQTL_all")

#result_df_total <- rbind(result_df_model1,result_df_model0)
result_df_total <- rbind(result_df_total,result_df_model3)
#result_df_total <- result_df_model2

#coef(fit_2, s=lambda.best) -> beta
#beta <- as.numeric(beta)
#rownames(coef(fit_2, s=lambda.best)) -> name_coeff
extract.coef(result_df_model3) -> test
test <- data.frame(test)
test$t <- test$Value/test$SE
test$p <- (1-pnorm(abs(test$t)))*2


result_dir <- paste0("/gpfs/ysm/pi/zhao-data/wl382/cWAS/data/fit_expr/",chr,"/",gene,"/")
dir.create(result_dir,recursive=TRUE,showWarnings =FALSE)
result_file <- paste0(result_dir,tissue,"_",tissue_frac,"_fitted_expr_r2_bothM0M1_M0M1.txt")
write.table(result_df_total,file=result_file,col.names=T,row.names=F,sep="\t",quote=F)
save(test,file=gsub("txt","Rdata",result_file))


