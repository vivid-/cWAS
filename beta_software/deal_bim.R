#!/usr/bin/R

args <- commandArgs()
#chr <- as.numeric(args[6])

result_df <- data.frame()
for(i in 1:22){
chr <- i
loc_file <- paste0("/ysm-gpfs/pi/zhao/ml2376/createdb/snplist/snplist.chr",chr,".txt")
loc <- read.table(loc_file,header=F,stringsAsFactors=F)
bim <- read.table("/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/GTEX/v6/geno/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.bim",header=F,stringsAsFactors=F)
bim_1 <- bim[bim[,1]==chr,]

loc_1 <- loc[loc[,1]==paste0("chr",chr),]
rsid <- loc_1[match(bim_1[,4],loc_1[,2]),3]

bim_1[!is.na(rsid),2] <- rsid[!is.na(rsid)]
result_df <- rbind(result_df,bim_1)
}

write.table(result_df,file="/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/GTEX/v6/geno/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs_hg19.bim",col.names=F,row.names=F,sep="\t",quote=F)
