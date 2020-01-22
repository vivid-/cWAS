#!/usr/bin/R
library("snpStats")
#args <- commandArgs()
#expr_model_beta_file <- args[6]
args <- commandArgs(trailingOnly = TRUE)

out_f <- match.fun(strsplit(grep('--out*', args, value = TRUE), split = '=')[[1]][[2]])
expr_beta_f <- match.fun(strsplit(grep('--expr-beta*', args, value = TRUE), split = '=')[[1]][[2]])
S_f <- as.numeric(strsplit(grep('--signature-matrix*', args, value = TRUE), split = '=')[[1]][[2]])
Z_f <- match.fun(strsplit(grep('--gwas-beta*', args, value = TRUE), split = '=')[[1]][[2]])
ref_geno_f <- match.fun(strsplit(grep('--ref-geno*', args, value = TRUE), split = '=')[[1]][[2]])
ref_frac_f <- match.fun(strsplit(grep('--ref-cellFraction*', args, value = TRUE), split = '=')[[1]][[2]])
tol <- match.fun(strsplit(grep('--tol-inv*', args, value = TRUE), split = '=')[[1]][[2]])
# quality control like removing ambiguous SNP
qc_snps <- function(snp_df){
#input format: SNP REF ALT
  ambi_list <- c("AT","CG","TA","CG")
  # connecting ref and alt
  ref_alt <- paste0(snp_df$REF,snp_df$ALT)
  # identify row indices for ambiguous SNPs
  ambi_ind <- which(ref_alt %in% ambi_list)
  
  return(ambi_ind)
}


# match SNPs in expr_beta, ref_geno and gwas_beta
# like matching the effect alleles in expr_beta and gwas_beta
match_snps <- function(expr_beta,ref_snp,ref_geno,gwas_beta){
  # find shared part among the three datasets
  expr_beta_identifier <- paste0(expr_beta[,1],"_",expr_beta[,2],"_",expr_beta[,3])
  expr_beta_identifier <- c(expr_beta_identifier,paste0(expr_beta[,1],"_",expr_beta[,3],"_",expr_beta[,2]))
  
  gwas_beta_identifier <- paste0(gwas_beta[,1],"_",gwas_beta[,2],"_",gwas_beta[,3])
  gwas_beta_identifier <- c(gwas_beta_identifier,paste0(gwas_beta[,1],"_",gwas_beta[,3],"_",gwas_beta[,2]))

  ref_snp_identifier <- paste0(ref_snp[,2],"_",ref_snp[,5],"_",ref_snp[,6])
  ref_snp_identifier <- c(ref_snp_identifier, paste0(ref_snp[,2],"_",ref_snp[,6],"_",ref_snp[,5]))

  inter_1 <- intersect(expr_beta_identifier,gwas_beta_identifier)
  inter_2 <- c(intersect(inter_1,ref_snp_identifier))
  # add ref and alt switched identifier
  inters <- switch_RA(inter_2)
  
  ref_snp_identifier <- paste0(ref_snp[,2],"_",ref_snp[,5],"_",ref_snp[,6])
  ref_snp_1 <- ref_snp[which(ref_snp_identifier %in% inters),]
  ref_geno_1 <- ref_geno[,which(ref_snp_identifier %in% inters)]
  
  # deal with the gwas_beta
  gwas_beta_1 <- match_summary(gwas_beta,ref_snp_1,4)
  expr_beta_1 <- match_summary(expr_beta,ref_snp_1,4:(dim(expr_beta)[2]))
  
  
  result<- list(ref_snp = ref_snp_1,ref_geno=ref_geno_1, gwas_beta= gwas_beta_1, expr_beta = expr_beta_1)
  return(result)
}

# deal with summary stats and switch REF and ALT if necessarily
match_summary <- function(gwas,ref_snp,beta_ind){
  ref_snp_identifier <- paste0(ref_snp[,2],"_",ref_snp[,5],"_",ref_snp[,6])
  ref_snp_identifier_reverse <- paste0(ref_snp[,2],"_",ref_snp[,6],"_",ref_snp[,5])
  gwas_identifier <- paste0(gwas[,1],"_",gwas[,2],"_",gwas[,3])
  result_df <- data.frame(matrix(NA,ncol=ncol(gwas),nrow=nrow(ref_snp)))
  ind <- match(ref_snp_identifier,gwas_identifier)
  result_df[!is.na(ind),] <- gwas[ind[!is.na(ind)],]
  
  ind_reverse <- match(ref_snp_identifier_reverse,gwas_identifier) 
  tmp <- gwas[ind_reverse[!is.na(ind_reverse)],]
  tmp[,2] <- gwas[ind_reverse[!is.na(ind_reverse)],3]
  tmp[,3] <- gwas[ind_reverse[!is.na(ind_reverse)],2]
  tmp[,beta_ind] <- -tmp[,beta_ind]
  result_df[is.na(ind_reverse)] <- tmp 
 
  return(result_df)
}

# switch the REF and ALT
switch_RA <- function(identifiers){
  result <- identifiers
  for(iden in identifiers){
    tmp <- strsplit(iden,"_")[[1]]
    tmp_1 <- c(tmp[1],"_",tmp[3],"_",tmp[2])
    if(tmp_1 %in% identifiers){
      next()
    }else{
      result <- c(result, tmp_1)
    }
  }

  return(result)
}


# calculate A
# tol: tolerance for singular values larger than tol are considered non-zero
calc_A <- function(signature,tol){
  s <- signature[,-1]
  tmp <- s %*% t(s)
  h <- try(chol2inv(chol(tmp)))
  if(inherits(h,'try-error')){
   h <- pseudoinverse(tmp, tol)
  }
  
  return(h)
}

# get the variance of cell fraction for each cell type in reference database
# calculate the vector for each cell  M_c=expr_beta %*% Signature (G\times C) %*% A_c 
calc_Mc <- function(expr_beta,signature,A){
  result_df <- data.frame
  for(i in 1:ncol(A)){
    a <- expr_beta[,-c(1:3)]
    b <- signature[,-c(1)]
    tmp <- a%*% t(b) %*% A[,i]
    result_df <- rbind(result_df,tmp)
  }

  return(result_df)
}


# read expr_beta_f # only for matched genes in the signature matrix
# format: SNP.ID REF ALT GENE1 GENE2 GENE3...
expr_beta <- read.table(expr_beta_f,header=T,stringsAsFactors=F,sep="\t")
# read signatrue matrix in the corresponding tissue
# format: CELL.type GENE1 GENE2 GENE3 ...
signature <- read.table(S_f,header=T,stringsAsFactors=F,sep="\t")

# match the genes in sigature and expr_beta
inter_genes <- intersect(colnames(expr_beta)[-1], colnames(signature)[-1])
expr_beta_1 <- expr_beta[,c(1,match(inter_genes,colnames(expr_beta)))]
signature_1 <- signature[,c(1,match(inter_genes,colnames(signature)))]

# read the gwas betas
# format: SNP.ID REF ALT Z N
gwas_beta <- read.table(Z_f,header=T,stringsAsFactors=F,sep="\t")
# read the reference genotype
# format: plink bim/bed/fam file
ref_geno <- read.plink(paste0(ref_geno_f,".bed"), paste0(ref_geno_f,".bim"),paste0(ref_geno_f,".fam"), na.strings = ("-9"))
# read the reference cell fraction
# format: Sample.ID Cell1 Cell2 ...
ref_frac <- read.table(ref_frac_f,header=T,stringsAsFactors=F,sep="\t")


##### remove ambiguous SNPs #####
# for expr_beta
ind_beta <- qc_snps(expr_beta_1[,c(1,2,3)])
expr_beta_2 <- expr_beta_1[ind_beta,]
# for gwas_beta
ind_gwas <- qc_snps(gwas_beta[,c(1,2,3)])
gwas_beta_2 <- gwas_beta[ind_gwas,]
# for ref_geno
ind_ref <- qc_snps(ref_geno$map[,c(2,5,6)])
qc_snps_1 <- qc_snps[ind_ref,]
ref_geno_1 <- ref_geno$genotype[,ind_ref]

# match the SNPs in expr_beta_1, gwas_beta and ref_geno
result_dfs <- match_snps(expr_beta_2,qc_snps_1,ref_geno_1,gwas_beta_2)
qc_snps_3 <- result_dfs[['ref_snp']]
ref_geno_3 <- result_dfs[["ref_geno"]]
expr_beta_3 <- result_dfs[["expr_beta"]]
gwas_beta_3 <- result_dfs[["gwas_beta"]]

# calculate the variance of genotypes
geno_se <- apply(ref_geno_3, 2, se)
frac_se <- apply(ref_frac[,-1],2,se)

A <- calc_A(signature_1,as.numeric(tol))
M <- calc_Mc(expr_beta_3,signature_1,A)
# for each cell type, get the test statistic
assoc_df <- data.frame()
for(i in 1:ncol(A)){
  tmp_1 <- gwas_beta_3[,"Z"] * M[i,]
  tmp_2 <- tmp_1 * geno_se
  z <- sum(tmp_2)/frac_se[i]
  p_tmp <- (1- pnorm(abs(z)))*2

  tmp_assoc<- data.frame(cell = signature_1[i,1],z = z, p=p_tmp)
  assoc_df <- rbind(assoc_df,tmp_assoc)
}

write.table(assoc_df,file=out_f,col.names=T,row.names=F,sep="\t",quote=F)
