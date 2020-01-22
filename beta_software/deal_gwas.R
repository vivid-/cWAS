#!/usr/bin/R

gwas_file <- "/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/gwas/longevity/JorisDeelen_2019_longevity_99th_percentile_hg19.txt"
gwas <- read.table(gwas_file,header=T,stringsAsFactors=F,sep="\t")

loc2id <- function(gwas_df,snp_df){
  id <- snp_df[match(gwas_df[,3],snp_df[,2]),3]
  return(id)
}

result_df <- data.frame()
for(i in 1:22){
  gwas_tmp <- gwas[gwas[,2]==i,]
  loc_file <- paste0("/ysm-gpfs/pi/zhao/ml2376/createdb/snplist/snplist.chr",i,".txt")
  loc <- read.table(loc_file,header=F,stringsAsFactors=F)
  
  gwas_tmp_id <- loc2id(gwas_tmp,loc)
  gwas_tmp$id <- gwas_tmp_id
  result_df <- rbind(result_df,gwas_tmp)
}

rsID <- result_df$id
ref <- toupper(result_df$NEA)
alt <- toupper(result_df$EA)
z <- result_df$Beta/result_df$SE
n <- result_df$Effective_N

gwas_df <- data.frame(rsID = rsID, ref=ref, alt=alt, Z=z, N=n)
gwas_df <- na.omit(gwas_df)
write.table(gwas_df,file="/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/gwas/longevity/JorisDeelen_2019_longevity_99th_percentile_hg19_curated.txt",col.names=T,row.names=F,sep="\t",quote=F)
