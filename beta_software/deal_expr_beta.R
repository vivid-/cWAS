#!/usr/bin/R

args <- commandArgs()
chr <- as.numeric(args[6])

bim <- read.table("/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/GTEX/v6/geno/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.txt",header=F,stringsAsFactors=F,colClasses = c("numeric","character","numeric","numeric","character","character","character"))
#bim <- read.table("/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/GTEX/v6/geno/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.bim",header=F,stringsAsFactors=F,colClasses = c("numeric","character","numeric","numeric","character","character"))
# extract bim in the chromosome
bim_1 <- bim[bim[,1]==chr,]
#bim_1 <- bim
expr_dir <- paste0("/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/GTEX/v6/weights/Whole_Blood/chr",chr,"/")
files <- list.files(expr_dir,pattern="_weights\\.txt")

result_df <- data.frame(rsID=bim_1[,2],ref=bim_1[,5],alt=bim_1[,6])
for(file in files){
  #cat(file,"\n")
  gene <- gsub("_weights\\.txt","",file)
  tmp_file <- paste0(expr_dir,file)
  cat(tmp_file,"\n")
  tmp <- read.table(paste0(expr_dir,file),header=F,stringsAsFactors=F,colClasses = c("character","character","numeric","character","character"))
  weight <- rep(0,length(bim_1[,2]))
  # matched indicies
  #a <- match(tmp[,1],result_df$rsID)
  a <- match(tmp[,1],bim_1[,7])
  weight[a[!is.na(a)]] <- tmp[!is.na(a),3]
  result_df[a[!is.na(a)],]$ref <- tmp[!is.na(a),4]
  result_df[a[!is.na(a)],]$alt <- tmp[!is.na(a),5]  
  tmp_df <- data.frame(weight=weight)

  coln_original <- colnames(result_df)
  result_df <- cbind(result_df,tmp_df)
  colnames(result_df) <- c(coln_original,gene)
}

#tmp <- strsplit(result_df$rsID,"_")
#chr <- sapply(tmp,"[",1)
#loc <- sapply(tmp,"[",2)
#result_df$rsID <- paste0(chr,":",loc)
#result_df$rsID <- gsub("_b37","",result_df$rsID)

write.table(result_df,file=paste0("/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/GTEX/v6/weights/Whole_Blood/expr_beta_chr",chr,".txt"),col.names=T,row.names=F,sep="\t",quote=F)
