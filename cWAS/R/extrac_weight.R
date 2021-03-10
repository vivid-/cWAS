#!/usr/bin/R

#' function to extract imputation weights for signature genes
#'
#' @title extrac_weight
#'
#' @param signature_f file location for signature matrix
#' @param gene_df_f file location for the file including gene ensembl ids and gene ids
#' @param weight_f file location for imputation weights
#' @param r2_f file location for imputation R2
#'
#' @export

extrac_weight <- function(signature_f,gene_df_f,weight_f,r2_f){

  matched_weights <- list()
  # load signature matrix
  signature <- readRDS(file=signature_f)
  genes <- colnames(signature)

  # matrix for gene ids and gene symbols
  gene_df <- read.table(gene_df_f,stringsAsFactors=F)
  gene_df <- gene_df[match(genes,gene_df[,4]),]
  gene_df_id <- gsub("\\.[0-9]+","",gene_df[,7])

  # extract genes whose imputation r2 is above median level
  r2 <- read.table(r2_f,header=T,stringsAsFactors = F,sep="\t")
  r2_sig <- r2[which(r2$rsq > median(r2$rsq)&r2$pval< 0.05/nrow(r2)),]

  # read weights
  weight_all <- read.table(weight_f,header=T,stringsAsFactors=F,sep="\t",
                           colClasses=c("character","character","numeric","numeric","character",
                                        "character","numeric","character","character"))
  weight_all$ref[which(weight_all$ref=="TRUE")] <- "T"
  weight_all$alt[which(weight_all$alt=="TRUE")] <- "T"
  weight_all <- weight_all[which(weight_all$gene%in%r2_sig$gene),]
  weight_id <- gsub("\\.[0-9]+","",weight_all$gene)
  result_df <- gene_df[gene_df_id%in%weight_id,]

  # for each chromosome
  for(i in 1:22){
    total <- data.frame()
    tmp <- result_df[result_df[,1]==paste0("chr",i),]
    gs <- tmp[,7]
    gs <- na.omit(gs)
    if(length(gs)==0){
      next()
    }

    for(g in gs){
      g_1 <- gsub("\\.[0-9]+","",g)
      tmp <- weight_all[grep(g_1,weight_all$gene),c(2,7,5,6)]
      gene_symbol <- result_df[grep(gsub("\\.[0-9]+","",g),result_df[,7]),4]
      colnames(tmp) <- c("rsID",gene_symbol,"ref","alt")
      tmp_1 <- data.frame(rsID=tmp$rsID,ref=tmp$ref,alt=tmp$alt,gene_symbol=tmp[,2])
      colnames(tmp_1)[4] <- gene_symbol
      if(length(total)==0){
        total <- tmp_1
      }else{
        snp_df <- tmp_1[,c("rsID","ref","alt")]
        gene_existing <- colnames(total)[!colnames(total) %in% c("rsID","ref","alt")]
        mm <- matrix(0,dim(snp_df)[1], length(gene_existing))
        snp_df <- cbind(snp_df,data.frame(mm))
        colnames(snp_df) <- c("rsID","ref","alt",gene_existing)
        gene_adding <- colnames(tmp_1)[!colnames(tmp_1) %in% c("rsID","ref","alt")]
        ur <- matrix(0,dim(total)[1], length(gene_adding))

        snp_df <- cbind(snp_df,tmp_1[,!colnames(tmp_1) %in% c("rsID","ref","alt")])
        colnames(snp_df) <- c("rsID","ref","alt",gene_existing,gene_adding)
        coln <- colnames(total)

        total <- cbind(total,data.frame(ur))
        colnames(total) <- c(coln,gene_adding)
        total <- rbind(total,snp_df)
      }
    }
    if(length(total)==0){
      next()
    }
    final <- data.frame()
    identifers <- unique(paste0(total[,1],"_",total[,2],"_",total[,3]))
    all <- paste0(total[,1],"_",total[,2],"_",total[,3])
    for(idens in identifers){
      tmp <- total[which(all==idens),]
      tmp_1 <- tmp[1,]
      if(nrow(tmp)!=1){
        tmp_1[,-c(1:3)] <- colSums(tmp[,-c(1:3)])
      }
      final <- rbind(final,tmp_1)
    }

    matched_weights[[i]] <- final
  }
  return(matched_weights)
}
