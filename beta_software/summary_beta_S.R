#!/usr/bin/R

args <- commandArgs()
wkd_dir <- args[6]
result_df <- data.frame()
for(i in 1:22){
  tmp_file <- paste0(wkd_dir,"/expr_beta_S_chr",i,".txt")
  if(!file.exists(tmp_file)){
     next()
  }
  cat(tmp_file,"\n")
  tmp <- read.table(tmp_file,header=T,stringsAsFactors=F,sep="\t")
  result_df <- rbind(result_df,tmp)
}
write.table(result_df,file=paste0(wkd_dir,"/expr_beta_S_all.txt"),col.names=T,row.names=F,sep="\t",quote=F)
