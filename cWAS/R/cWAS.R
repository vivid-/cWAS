#' cWAS is a tool to identify cell types whose
#' regulated genetic components are disease-associated
#' library("genio")
#' library("MASS)
#'
#' @title cWAS
#'
#' @author Wei Liu
#'
#' @rdname cWAS
#'
#' @aliases cWAS
#'
#' @export
#'
#' @param gene_df_f file location for the file including gene ensembl ids and gene ids
#' @param weight_f file location for imputation weights
#' @param r2_f file location for imputation R2
#' @param signature_rds_f  file location for signature matrix rds file
#' @param gwas_f file location for gwas summary stats
#' @param ref_geno_f file location for reference genome
#' @param output_f output file prefix location
#'



cWAS <- function(gene_df_f,weight_f,r2_f,
                 signature_f,
                 gwas_f,ref_geno_f,output_f) {
  library(data.table)
  # match the signature gene and imputation weights to get the matched file
  matched_weights = extrac_weight(signature_f,gene_df_f,weight_f,r2_f)
  signature <- readRDS(file=signature_f)
  # for each chromosome
  # match SNPs
  matched_lists = list()
  # read the gwas betas
  # format: SNP.ID REF ALT Z N
  Z_f = gwas_f
  cat("loading gwas file......\n")
  if(length(grep("\\.gz$",Z_f,value=T))!=0){
    gwas_beta <- data.frame(fread(paste0("zcat ",Z_f),header=T,stringsAsFactors=F))
  }else{
    gwas_beta <- data.frame(fread(Z_f,header=T,stringsAsFactors=F))
  }

  # SNP	A2	A1	Z	N
  if("Z" %in% c(colnames(gwas_beta))){
    gwas_beta <- gwas_beta[,c("SNP","A2","A1","Z","N")]
  }else{
    gwas_beta$Z <- sign(as.numeric(gwas_beta$BETA))*(1-qnorm(as.numeric(gwas_beta$P)/2))
    gwas_beta <- gwas_beta[,c("SNP","A2","A1","Z","N")]
  }
  cat("gwas file loaded......\n")
  # read the reference genotype
  # format: plink bim/bed/fam file
  ref_geno <- read.table(paste0(ref_geno_f,".bim"),header=F,stringsAsFactors=F,sep="\t")

  cat("SNP and gene matching......\n")
  for(i in 1:22){
    cat("starting chr",i,"\n")
    expr_beta_tmp = matched_weights[[i]]
    if(is.null(expr_beta_tmp)){
      matched_lists[[i]] = NULL
      next()
    }
    tmp <- match_gwas_expr_ref(expr_beta_tmp,signature,gwas_beta,ref_geno)
    matched_lists[[i]] = tmp
    cat("finishing chr",i,"\n")
  }
  cat("SNP and gene matched......\n")
  # for each chromosome
  # calcuate beta_S

  cat("summarizing results......\n")
  beta_S_list <- list()
  for(i in 1:22){
    matched_tmp = matched_lists[[i]]
    if(is.null(matched_tmp)){
      beta_S_list[[i]] = NULL
      next()
    }
    expr_beta_tmp = matched_tmp[["expr_beta"]]
    sig_tmp = matched_tmp[["signature"]]
    tmp = get_expr_beta_s(expr_beta_tmp,sig_tmp)
    beta_S_list[[i]] = tmp
  }

  # concatenate all by-chr results
  expr_beta_all = summary_expr_beta(matched_lists)

  # concatenate beta_S results
  beta_S_all = data.frame()
  for(i in 1:22){
    tmp = beta_S_list[[i]]
    if(is.null(tmp)){
      next()
    }
    beta_S_all <- rbind(beta_S_all,tmp)
  }

  # concatenate gwas and ref_snp results
  matched_all <- summary_others(matched_lists)
  S_all = matched_all[["signa"]]
  gwas_all <- matched_all[["gwas"]]
  bim_all <- matched_all[['bim']]

  cat("printing results......\n")
  # extract the corresponding reference genotype
  write.table(bim_all,file=paste0(output_f,".bim"),col.names=T,row.names=F,sep="\t",quote=F)
  write.table(gwas_all,file=paste0(output_f,".gwas"),col.names=T,row.names=F,sep="\t",quote=F)
  write.table(S_all,file=paste0(output_f,".signature"),col.names=T,row.names=F,sep="\t",quote=F)
  write.table(beta_S_all,file=paste0(output_f,".betaS"),col.names=T,row.names=F,sep="\t",quote=F)
  write.table(expr_beta_all,file=paste0(output_f,".beta"),col.names=T,row.names=F,sep="\t",quote=F)
}


#'function
#' @title cWAS_test
#' @param ref_geno_f file location for after-match reference genome
#' @param S_f file location for after-match signature matrix
#' @param Z_f file location for after-match gwas
#' @param betaS_f file location for after-match betaS
#' @param out_f  file prefix for file location
#'
#' @export
#'
#'library("genio")
#'library("MASS")
#'
cWAS_test <- function(ref_geno_f,S_f,Z_f,
                      betaS_f,out_f){
  library(genio)
  # assuming all data (expr_beta, ref_geno, Z and S) are all matched
  # read signatrue matrix in the corresponding tissue
  # format: CELL.type GENE1 GENE2 GENE3 ...
  signature <- read.table(S_f,header=T,stringsAsFactors=F,sep="\t")
  # determine the number of genes
  if(ncol(signature)<50){
    stop("Number of signature genes fewer than 50")
  }

  # read the gwas betas
  # format: SNP.ID REF ALT Z N
  gwas_beta <- read.table(Z_f,header=T,stringsAsFactors=F,sep="\t")
  inds <- which(!is.na(gwas_beta$Z))

  ref_geno <- read_plink(ref_geno_f)
  bim <- data.frame(ref_geno$bim)
  bim_all_ind <- paste0(gwas_beta[,1],"_",gwas_beta[,2],"_",gwas_beta[,3])
  bim_all_ind <- bim_all_ind[inds]
  bim_ind <- paste0(bim$id,"_",bim$ref,"_",bim$alt)
  ind_match <- match(bim_all_ind,bim_ind)
  bim <- bim[ind_match,]
  # read the reference cell fraction
  # format: Sample.ID Cell1 Cell2 ...
  # ref_frac <- read.table(ref_frac_f,header=T,stringsAsFactors=F,sep="\t")
  # calculate the variance of genotypes
  # deal with geno
  geno <- ref_geno$X
  geno <- geno[ind_match,]
  geno_se <- apply(geno, 1, function(X) sd(na.omit(X)))
  geno_mean <- apply(geno, 1, function(X) mean(na.omit(X)))
  # extract matched SNPs
  inds_na <- which(is.na(geno),arr.ind=T)
  geno[inds_na] <- geno_mean[inds_na[,1]]
#  geno <- geno[ind_match,]
  gwas_beta <- gwas_beta[inds,]
  geno <- geno[inds,]

  # read betas
  betaS <- read.table(betaS_f,header=T,stringsAsFactors=F,sep="\t")
  betaS <- read.table(betaS_f,header=T,stringsAsFactors=F,sep="\t",
                      colClasses = c(rep("numeric",ncol(betaS)-3),
                                     "character", "character", "character"))

  betaS$ref[which(betaS$ref=="TRUE")] <- "T"
  betaS$alt[which(betaS$alt=="TRUE")] <- "T"
  betaS <- betaS[inds,]
  betaS_iden <- paste0(betaS$rsID,"_",betaS$ref,"_",betaS$alt)
  #bim <- data.frame(ref_geno$bim)
  # extract matched SNPs
  #bim <- bim[,]
  bim <- bim[inds,]
  bim_iden <- paste0(bim$id,"_",bim$ref,"_",bim$alt)
  inter_iden <- intersect(betaS_iden,bim_iden)
  geno <- geno[match(inter_iden,bim_iden),]
  betaS <- betaS[match(inter_iden,betaS_iden),]
  geno_se <- apply(geno, 1, function(X) sd(na.omit(X)))
  geno_mean <- apply(geno, 1, function(X) mean(na.omit(X)))
  betaS <- betaS[,-c(ncol(betaS)-2,ncol(betaS)-1,ncol(betaS))]
  A <- calc_A(signature,as.numeric(tol))

  f <- t(geno) %*% as.matrix(betaS) %*% A
  colnames(f) <- signature[,1]
  frac_se <- apply(f,2,sd)
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

  #bim_inter <- bim[match(inter_iden,bim_iden),]
  #snp_cell <- data.frame(rsid=bim_inter$id,ref=bim_inter$ref,alt=bim_inter$alt)
  assoc_df_sig <- assoc_df[assoc_df$p <0.05/nrow(assoc_df),]
  ind_sig <- which(assoc_df$p <0.05/nrow(assoc_df))
  assoc_df <- assoc_df[order(assoc_df$p),]

  write.table(assoc_df,file=out_f,col.names=T,row.names=F,sep="\t",quote=F)
}


#'function to concatenate reference genome and gwas
#'
#'@param matched_list the output from match_gwas_expr_ref
#'
#'@export
#'
summary_others <- function(matched_list){
  gwas <- data.frame()
  #expr_beta <- data.frame()
  signa <- data.frame()
  bim <- data.frame()
  for(i in 1:22){
    if(is.null(matched_list[[i]])){
      next()
    }
    gwas <- rbind(gwas,matched_list[[i]][["gwas_beta"]])
    bim <- rbind(bim,matched_list[[i]][['ref_snp']])

    sig_tmp <- matched_list[[i]][["signature"]]
    if(nrow(signa)==0){
      signa = sig_tmp
    }else{
      coln_o = colnames(signa)
      signa = cbind(signa,sig_tmp[,-1])
      colnames(signa) = c(coln_o,colnames(sig_tmp)[-1])
    }
  }

  result_list <- list()
  result_list[['gwas']] = gwas
  result_list[['bim']] = bim
  result_list[['signa']] = signa
  return(result_list)
}


#'function to concatenate expr_beta
#'
#'@param matched_list the output from match_gwas_expr_ref
#'
#'@export
#'

summary_expr_beta <- function(matched_list){
  result_df <- data.frame()
  for(i in 1:22){
    chr = i
    tmp = matched_list[[i]]
    if(is.null(tmp)){
      next
    }

    tmp = tmp[["expr_beta"]]
    if(nrow(result_df)==0){
      result_df <- tmp
    }else{
      # for existing gene add snps
      snp_df <- tmp[,c("rsID","ref","alt")]
      gene_existing <- colnames(result_df)[!colnames(result_df) %in% c("rsID","ref","alt")]
      mm <- matrix(0,dim(snp_df)[1], length(gene_existing))
      snp_df <- cbind(snp_df,data.frame(mm))
      colnames(snp_df) <- c("rsID","ref","alt",gene_existing)

      # zeros in the right upper corner
      gene_adding <- colnames(tmp)[!colnames(tmp) %in% c("rsID","ref","alt")]
      ur <- matrix(0,dim(result_df)[1], length(gene_adding))
      snp_df <- cbind(snp_df,tmp[,!colnames(tmp) %in% c("rsID","ref","alt")])
      colnames(snp_df) <- c("rsID","ref","alt",gene_existing,gene_adding)

      # add genes
      coln <- colnames(result_df)
      result_df <- cbind(result_df,data.frame(ur))
      colnames(result_df) <- c(coln,gene_adding)
      result_df <- rbind(result_df,snp_df)
    }
  }

  return (result_df)
}
