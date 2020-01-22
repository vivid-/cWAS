#!/usr/bin/R
args <- commandArgs()
#tmp_dir <- "/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/assoc/"
tmp_dir <- args[6]
gwas <- data.frame()
expr_beta <- data.frame()
signature <- data.frame()
bim <- data.frame()
for(i in 1:22){
  cat(i,"\n")
  gwas_tmp <- try(read.table(paste0(tmp_dir,"gwas_chr",i,".txt"),header=T,stringsAsFactors=F,sep="\t"))
  if(inherits(gwas_tmp,"try-error")){
    next()
  }
  
  gwas <- rbind(gwas,gwas_tmp)

#  expr_tmp <- read.table(paste0(tmp_dir,"WB_longevity_expr_beta_v6_chr",i,".txt"),header=T,stringsAsFactors=F,sep="\t")
#  if(length(signature)==0){
#    expr_beta <- expr_tmp
#  }else{
#    coln_o <- colnames(expr_beta)
i#    expr_beta <- cbind(expr_beta,expr_tmp[,-c(1:3)])
#
#    colnames(expr_beta) <- c(coln_o,colnames(expr_tmp)[-c(1:3)])   
#  } 
  #expr_beta <- rbind(expr_beta,expr_tmp)
  
  bim_tmp <- read.table(paste0(tmp_dir,"gtex_v8_matched_chr",i,".bim"),header=F,stringsAsFactors=F,sep="\t")
  bim <- rbind(bim,bim_tmp)

  sig_tmp <- read.table(paste0(tmp_dir,"LM22_chr",i,".txt"),header=T,stringsAsFactors=F,sep="\t")
  if(length(signature)==0){
    signature <- sig_tmp
  }else{
    coln_o <- colnames(signature)
    signature <- cbind(signature,sig_tmp[,-1])
    colnames(signature) <- c(coln_o,colnames(sig_tmp)[-1])
  }
}

write.table(gwas,file=paste0(tmp_dir,"gwas_all.txt"),col.names=T,row.names=F,sep="\t",quote=F)
#write.table(expr_beta,file=paste0(tmp_dir,"gtex_v6_matched_all.txt"),col.names=T,row.names=F,sep="\t",quote=F)
write.table(signature,file=paste0(tmp_dir,"LM22_new_all.txt"),col.names=T,row.names=F,sep="\t",quote=F)
write.table(bim,file=paste0(tmp_dir,"gtex_v8_all.bim"),col.names=F,row.names=F,sep="\t",quote=F)
