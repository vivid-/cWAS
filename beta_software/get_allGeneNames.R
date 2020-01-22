#!/usr/bin/R

col_names <- c()
for(i in 1:22){
  cat(i,"\n")
  chr <- i
  file=paste0("/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/GTEX/v6/weights/Whole_Blood/expr_beta_chr",chr,".txt")
  tmp <- read.table(file,header=T,stringsAsFactors=F,sep="\t")
  
  if(i ==1){
    col_names <- colnames(tmp)
  }else{
    # for existing gene add snps
    #snp_df <- tmp[,c("rsID","ref","alt")]
    #gene_existing <- colnames(result_df)[!colnames(result_df) %in% c("rsID","ref","alt")]
    #mm <- matrix(0,dim(snp_df)[1], length(gene_existing))
    #snp_df <- cbind(snp_df,data.frame(mm))
    #colnames(snp_df) <- c("rsID","ref","alt",gene_existing)
    
    # append rows
    #result_df<- rbind(result_df,snp_df)
    
    # zeros in the right upper corner
    gene_adding <- colnames(tmp)[!colnames(tmp) %in% c("rsID","ref","alt")]
    #ur <- matrix(0,dim(result_df)[1], length(gene_adding))
    #add_df <- rbind(data.frame(ur),tmp[,!colnames(tmp) %in% c("rsID","ref","alt")])
    
    #snp_df <- cbind(snp_df,tmp[,!colnames(tmp) %in% c("rsID","ref","alt")])
   
    # add rows first
    ##result_df<- rbind(result_df,snp_df)
    
    # add genes
    #coln <- colnames(result_df)
    #result_df <- cbind(result_df,data.frame(ur))
    #colnames(result_df) <- c(coln,gene_adding)
    #result_df <- rbind(result_df,snp_df)
    col_names <- c(col_names,gene_adding)
  }
}

write.table(col_names,file="/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/GTEX/v6/weights/Whole_Blood/expr_beta_colnames.txt",col.names=F,row.names=F,sep="\t",quote=F)

