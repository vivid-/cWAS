#!/usr/bin/R
args <- commandArgs()
wrk_dir <- args[6]
result_df <- data.frame()
for(i in 1:22){
  chr <- i
  file=paste0(wrk_dir,"/expr_beta_v8_chr",chr,".txt")
  if(!file.exists(file)){
    next()
  }
  tmp <- read.table(file,header=T,stringsAsFactors=F,sep="\t")
  
  if(i ==1){
    result_df <- tmp
  }else{
    # for existing gene add snps
    snp_df <- tmp[,c("rsID","ref","alt")]
    gene_existing <- colnames(result_df)[!colnames(result_df) %in% c("rsID","ref","alt")]
    mm <- matrix(0,dim(snp_df)[1], length(gene_existing))
    snp_df <- cbind(snp_df,data.frame(mm))
    colnames(snp_df) <- c("rsID","ref","alt",gene_existing)
    
    # append rows
    #result_df<- rbind(result_df,snp_df)
    
    # zeros in the right upper corner
    gene_adding <- colnames(tmp)[!colnames(tmp) %in% c("rsID","ref","alt")]
    ur <- matrix(0,dim(result_df)[1], length(gene_adding))
    #add_df <- rbind(data.frame(ur),tmp[,!colnames(tmp) %in% c("rsID","ref","alt")])
    
    snp_df <- cbind(snp_df,tmp[,!colnames(tmp) %in% c("rsID","ref","alt")])
    colnames(snp_df) <- c("rsID","ref","alt",gene_existing,gene_adding)
    # add rows first
    #result_df<- rbind(result_df,snp_df)
    
    # add genes
    coln <- colnames(result_df)
    result_df <- cbind(result_df,data.frame(ur))
    colnames(result_df) <- c(coln,gene_adding)
    result_df <- rbind(result_df,snp_df)
   
  }
}
write.table(result_df,file=paste0(wrk_dir,"/expr_beta.txt"),col.names=T,row.names=F,sep="\t",quote=F)
