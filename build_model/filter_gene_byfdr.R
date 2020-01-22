#!/usr/bin/R
args <- commandArgs()
tissue <- args[6]
#df <- read.table("/ysm-gpfs/pi/zhao-data/wl382/cWAS/scripts/git_up/build_model/loc_netid.txt",header=T,stringsAsFactors=F,sep="\t")
#cov_fl = list.files("/ysm-gpfs/pi/zhao-data/zy92/GTEx_V8/GTEx_Analysis_v8_eQTL_covariates/")
#idf = matrix(unlist(strsplit(cov_fl, "\\.")), ncol = 4, byrow = T)[,1]
##print(paste0("INFO: idf ", idf))
#tissues = as.character(sapply(cov_fl, function(x) unlist(strsplit(x, '\\.'))[1]))
netids <- c("wl382")
result_df <- data.frame()
for(chr in 1:22){
#  netids <- df[df[,1]==chr,2]
  for(netid in netids){
  genes <- list.files(paste0("/ysm-gpfs/scratch60/",netid,"/GTEx/v8/models/chr",chr,"/")) 
#  for(t in tissues){
  for(gene in genes){
  	paste0("/ysm-gpfs/scratch60/",netid,"/GTEx/v8/models/chr",chr,"/",gene,"/",tissue,"_rsq.txt") -> tmp
        dat <- try(read.table(tmp,header=T,stringsAsFactors=F))
        if(!inherits(dat,"try-error")){
           dat$chr <- chr
           result_df <- rbind(result_df,dat)
        }
#  } 
  }
  }
}

result_df <- unique(result_df)
write.table(result_df,file=paste0("/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/model/r2/",tissue,"_all_summary.txt"),col.names=T,row.names=F,sep="\t",quote=F)
