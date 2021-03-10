#!/usr/bin/R
#library(snpStats)
#library(data.table)

#' function to match the snps existing in gwas summary stats,
#' imputation weights and reference genome files
#'
#' @title match_gwas_expr_ref
#'
#' @param expr_beta imputation weights for signature genes
#' @param S_f file location for signature matrix
#' @param Z_f file location of gwas summary stats
#' @param ref_geno_f file location of reference genome
#'
#' @export

match_gwas_expr_ref <- function(expr_beta,signature,gwas_beta,ref_geno){
  # read signatrue matrix in the corresponding tissue
  # format: CELL.type GENE1 GENE2 GENE3 ...
  #signature <- read.table(S_f,header=T,stringsAsFactors=F,sep="\t")

  # match the genes in sigature and expr_beta
  inter_genes <- intersect(colnames(expr_beta)[-c(1:3)], colnames(signature))
  expr_beta_1 <- expr_beta[,c(1:3,match(inter_genes,colnames(expr_beta)))]

  signature_1 <- signature[,c(match(inter_genes,colnames(signature)))]
  a <- colnames(signature_1)
  signature_1 <- cbind(data.frame(cell_type=rownames(signature)),signature_1)
  colnames(signature_1)[-1] = inter_genes
  if(length(inter_genes)==0){
    cat("stop!")
    next()
  }


  expr_beta_2 <- expr_beta_1
  gwas_beta_2 <- gwas_beta
  qc_snps_1 <- ref_geno
  result_dfs <- match_snps(expr_beta_2,qc_snps_1,gwas_beta_2)
  if(!is.null(result_dfs)){
    qc_snps_3 <- result_dfs[['ref_snp']]
    expr_beta_3 <- result_dfs[["expr_beta"]]
    gwas_beta_3 <- result_dfs[["gwas_beta"]]
    result_dfs[["signature"]] = signature_1
  }

  return(result_dfs)
}


# quality control like removing ambiguous SNP
qc_snps <- function(snp_df){
  #input format: SNP REF ALT
  ambi_list <- c("AT","CG","TA","CG")
  # connecting ref and alt
  ref_alt <- paste0(snp_df[,2],snp_df[,3])
  # identify row indices for ambiguous SNPs
  ambi_ind <- which(ref_alt %in% ambi_list)

  return(ambi_ind)
}

# match SNPs in expr_beta, ref_geno and gwas_beta
# like matching the effect alleles in expr_beta and gwas_beta
match_snps <- function(expr_beta,ref_snp,gwas_beta){
  # find shared part among the three datasets
  expr_beta_identifier <- paste0(expr_beta[,1],"_",expr_beta[,2],"_",expr_beta[,3])
  expr_beta_identifier <- c(expr_beta_identifier,paste0(expr_beta[,1],"_",expr_beta[,3],"_",expr_beta[,2]))

  gwas_beta_identifier <- paste0(gwas_beta[,1],"_",gwas_beta[,2],"_",gwas_beta[,3])
  gwas_beta_identifier <- c(gwas_beta_identifier,paste0(gwas_beta[,1],"_",gwas_beta[,3],"_",gwas_beta[,2]))

  ref_snp_identifier <- paste0(ref_snp[,2],"_",ref_snp[,5],"_",ref_snp[,6])
  ref_snp_identifier <- c(ref_snp_identifier, paste0(ref_snp[,2],"_",ref_snp[,6],"_",ref_snp[,5]))

  inter_1 <- intersect(expr_beta_identifier,gwas_beta_identifier)
  inter_2 <- c(intersect(inter_1,ref_snp_identifier))
  if(length(inter_2)==0){
    return (NULL)
  }
  # add ref and alt switched identifier
  inters <- switch_RA(inter_2)

  ref_snp_identifier <- paste0(ref_snp[,2],"_",ref_snp[,5],"_",ref_snp[,6])
  ref_snp_1 <- ref_snp[which(ref_snp_identifier %in% inters),]
  #ref_geno_1 <- ref_geno[,which(ref_snp_identifier %in% inters)]
  colnames(ref_snp_1) <- colnames(ref_snp)
  # deal with the gwas_beta
  gwas_beta_1 <- match_summary(gwas_beta,ref_snp_1,4)
  expr_beta_1 <- match_summary(expr_beta,ref_snp_1,4:(dim(expr_beta)[2]))
  colnames(expr_beta_1) <- colnames(expr_beta)
  colnames(gwas_beta_1) <- colnames(gwas_beta)

  result<- list(ref_snp = ref_snp_1, gwas_beta= gwas_beta_1, expr_beta = expr_beta_1)
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
  result_df[is.na(ind),] <- tmp

  return(result_df)
}


# switch the REF and ALT
switch_RA <- function(identifiers){

  tmp <- strsplit(identifiers,"_")
  chr <- sapply(tmp,"[",1)
  ref <- sapply(tmp,"[",2)
  alt <- sapply(tmp,"[",3)
  tmp_1 <- paste0(chr,"_",alt,"_",ref)
  result <- c(unique(c(tmp_1,identifiers)))

  return(result)
}


