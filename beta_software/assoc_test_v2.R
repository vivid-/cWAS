#!/usr/bin/R

# after matching SNPs, getting beta%*%S, calculate the Z scores for each cell type
library("genio")
#args <- commandArgs()
#expr_model_beta_file <- args[6]
args <- commandArgs(trailingOnly = TRUE)

out_f <- strsplit(grep('--out*', args, value = TRUE), split = '=')[[1]][[2]]
S_f <- strsplit(grep('--signature-matrix*', args, value = TRUE), split = '=')[[1]][[2]]
Z_f <- strsplit(grep('--gwas-beta*', args, value = TRUE), split = '=')[[1]][[2]]
ref_geno_f <- strsplit(grep('--ref-geno*', args, value = TRUE), split = '=')[[1]][[2]]
ref_frac_f <- strsplit(grep('--ref-cellFraction*', args, value = TRUE), split = '=')[[1]][[2]]
tol <- as.numeric(strsplit(grep('--tol-inv*', args, value = TRUE), split = '=')[[1]][[2]])
betaS_f <- strsplit(grep('--betaS-file*', args, value = TRUE), split = '=')[[1]][[2]]

# calculate A
# tol: tolerance for singular values larger than tol are considered non-zero
calc_A <- function(signature,tol){
  s <- as.matrix(signature[,-1])
  tmp <- s %*% t(s)
  h <- try(chol2inv(chol(tmp)))
  if(inherits(h,'try-error')){
   h <- pseudoinverse(tmp, tol)
  }
  
  return(h)
}

# get the variance of cell fraction for each cell type in reference database
# calculate the vector for each cell  M_c=expr_beta %*% Signature (G\times C) %*% A_c 
calc_Mc <- function(betaS,A){
  result_df <- data.frame()
  A <- as.matrix(A)
  betaS <- as.matrix(betaS)
  for(i in 1:ncol(A)){
    #a <- expr_beta[,-c(1:3)]
    #b <- signature[,-c(1)]
    tmp <- betaS %*% A[,i]
    result_df <- rbind(result_df,t(tmp))
  }

  return(result_df)
}


# assuming all data (expr_beta, ref_geno, Z and S) are all matched
# read signatrue matrix in the corresponding tissue
# format: CELL.type GENE1 GENE2 GENE3 ...
signature <- read.table(S_f,header=T,stringsAsFactors=F,sep="\t")
# read the gwas betas
# format: SNP.ID REF ALT Z N
gwas_beta <- read.table(Z_f,header=T,stringsAsFactors=F,sep="\t")
# read the reference genotype
# format: plink bim/bed/fam file
#ref_geno <- read.plink(paste0(ref_geno_f,".bed"), paste0(ref_geno_f,".bim"),paste0(ref_geno_f,".fam"), na.strings = ("-9"))
ref_geno <- read_plink(ref_geno_f)
# read the reference cell fraction
# format: Sample.ID Cell1 Cell2 ...
ref_frac <- read.table(ref_frac_f,header=T,stringsAsFactors=F,sep="\t")
# calculate the variance of genotypes
geno_se <- apply(ref_geno$X, 1, function(X) sd(na.omit(X)))
frac_se <- apply(ref_frac[,-ncol(ref_frac)],2,sd)
# read \hat\beta %*% S
# format: Cell1 Cell2 Cell3...rsID
betaS <- read.table(betaS_f,header=T,stringsAsFactors=F,sep="\t")
betaS <- betaS[,-c(ncol(betaS)-2,ncol(betaS)-1,ncol(betaS))]
A <- calc_A(signature,as.numeric(tol))
#M <- calc_Mc(betaS,A)
M <- A %*% t(betaS)

assoc_df <- data.frame()
for(i in 1:ncol(A)){
  tmp_1 <- gwas_beta[,"Z"] * M[i,]
  tmp_2 <- tmp_1 * geno_se
  z <- sum(tmp_2)/frac_se[i]
  p_tmp <- (1- pnorm(abs(z)))*2

  tmp_assoc<- data.frame(cell = signature[i,1],z = z, p=p_tmp)
  assoc_df <- rbind(assoc_df,tmp_assoc)
}
#assoc_df <- assoc_df[order(assoc_df$p),]
#write.table(assoc_df,file=out_f,col.names=T,row.names=F,sep="\t",quote=F)

# only for significant cell types
assoc_df_sig <- assoc_df[assoc_df$p <0.05/nrow(assoc_df),]
ind_sig <- which(assoc_df$p <0.05/nrow(assoc_df))
assoc_df <- assoc_df[order(assoc_df$p),]
write.table(assoc_df,file=out_f,col.names=T,row.names=F,sep="\t",quote=F)

#A_sig <- calc_A(signature[ind_sig,],as.numeric(tol))
#M_sig <- A_sig %*% t(betaS[,ind_sig])
A_sig <- calc_A(signature,as.numeric(tol))
M_sig <- A_sig %*% t(betaS)
#f <- as.matrix(ref_frac[,ind_sig])
f <- as.matrix(ref_frac[,-ncol(ref_frac)])
#ind_na <- which(is.na(f),arr.ind=T)
# remove columns with all na
colSums(f,na.rm=TRUE) -> colsum
ind_col <- which(colsum==0)
rowSums(f,na.rm=TRUE) -> rowsum
ind_row <- which(rowsum==0)
if(length(ind_col)!=0){
  f <- f[,-ind_col]
  A_sig <- A_sig[-ind_col,-ind_col]
  M_sig <- M_sig[-ind_col,]
}
if(length(ind_row)!=0){
  f <- f[-ind_row,]
}
tmp <- t(f) %*% f

U <- try(chol2inv(chol(tmp)))
if(inherits(U,'try-error')){
U <- pseudoinverse(tmp, tol)
}


H <- U%*%M_sig
assoc_sig <- data.frame()
for(i in 1:ncol(A_sig)){
  tmp_1 <- gwas_beta[,"Z"] * H[i,]
  tmp_2 <- tmp_1 * geno_se
  #z <- sum(tmp_2)/frac_se[i]
  z <- sum(tmp_2)/sqrt(U[i,i])
  p_tmp <- (1- pnorm(abs(z)))*2

  tmp_assoc<- data.frame(cell = signature[i,1],z = z, p=p_tmp)
  assoc_sig <- rbind(assoc_sig,tmp_assoc)
}

assoc_sig <- assoc_sig[order(assoc_sig$p),]
write.table(assoc_sig,file=paste0(out_f,".conditional"),col.names=T,row.names=F,sep="\t",quote=F)
