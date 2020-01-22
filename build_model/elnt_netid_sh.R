locs <- read.table("loc_netid.txt",header=T,stringsAsFactors=F)

sink("r2_allgenes.sh")
for(i in 1:22){
  chr <- paste0("chr",i)
  netids <- locs[locs[,1]==i,2]
  for(netid in netids){
    genes <- list.files(paste0("/ysm-gpfs/scratch60/",netid,"/GTEx/v8/models/",chr,"/"))
    for(gene in genes){
      result_dir <- paste0("/ysm-gpfs/scratch60/wl382/GTEx/v8/models/",chr,"/",gene)
      cat(paste0("module load R/3.6.1-foss-2018b;mkdir -p ",result_dir,";Rscript /ysm-gpfs/pi/zhao-data/wl382/cWAS/scripts/git_up/build_model/elnt_netid_rsq.R chr",i," ",gene," ",netid),"\n")
    }
  }
}
sink()
