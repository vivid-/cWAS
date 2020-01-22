#!/usr/bin/R

#for(i in 1:22){
  
#}
args <- commandArgs()
chr <- args[6]

tmp_dir <- paste0("/gpfs/ysm/pi/zhao-data/zd72/cWAS/GTEx/v6/pred/ukbb/Lung/",chr,"/")
files <- list.files(tmp_dir,pattern="expr.profile")

tmp_file <- paste0(tmp_dir,files[1])
tmp <- read.table(tmp_file,header=T,stringsAsFactors=F)
result_df <- data.frame(IID=tmp$IID)
for(file in files){
  tmp_file <- paste0(tmp_dir,file)
  tmp <- read.table(tmp_file,header=T,stringsAsFactors=F)
  gene <- gsub("_ukbb_lung_pred.expr\\.profile","",file)
  #tmp_df <- data.frame(num = tmp$SCORESUM)
  tmp_df <- data.frame(num = rep(0,length(tmp$IID)))
  colnames(tmp_df) <- gene
  tmp_df[match(tmp$IID,result_df$IID),1] <- tmp$SCORESUM 
  

  result_df <- cbind(result_df,tmp_df)
}

outputf <- args[7]
write.table(result_df,file=outputf,col.names=T,row.names=F,sep="\t",quote=F)
