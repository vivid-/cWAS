#!/usr/bin/R
args <- commandArgs()
chr <- args[6]
netid <- args[7]
cov_fl = list.files("/ysm-gpfs/pi/zhao-data/zy92/GTEx_V8/GTEx_Analysis_v8_eQTL_covariates/")
idf = matrix(unlist(strsplit(cov_fl, "\\.")), ncol = 4, byrow = T)[,1]
#print(paste0("INFO: idf ", idf))
tissues = as.character(sapply(cov_fl, function(x) unlist(strsplit(x, '\\.'))[1]))

outputf <- paste0("models_chr",chr,"_",netid,"_remain.sh")
sink(outputf)

#for(chr in 1:22){
  chr <- paste0("chr",chr)
  genes = list.files(paste0("/ysm-gpfs/project/wl382/GTEx_v8/genotype/cis_loc/",chr,"/"))
  for(gene in genes){
    for(tissue in tissues){
       result_dir <- paste0("/ysm-gpfs/scratch60/",netid,"/GTEx/v8/models/",chr,"/",gene)
       gene_i <- gsub("\\.[0-9]+","",gene)
       expr_dir <- paste0("/gpfs/loomis/scratch60/fas/radev/zy92/GTEX/adjusted_expr1/",chr,"/",gene_i,"/",tissue,".adj_expr")
       if(!file.exists(expr_dir)){
          next()
       }
       cis_snp_dir <- paste0("/ysm-gpfs/project/wl382/GTEx_v8/genotype/cis_loc/",chr,"/",gene)
       if(!file.exists(cis_snp_dir)){
          next()
       }
       if(file.exists(paste0("/ysm-gpfs/scratch60/wl382/GTEx/v8/models/",chr,"/",gene,"/",tissue,"_weight.txt"))){
          next()

       }
       sent <- paste0("mkdir -p ",result_dir,";module load R/3.6.1-foss-2018b; Rscript /ysm-gpfs/pi/zhao-data/wl382/cWAS/scripts/git_up/build_model/elnt_netid.R ",chr," ",gene," ",tissue," ",netid)
       cat(sent,"\n")
    }     
  
  }
#}
sink()
