#!/usr/bin/R
cov_fl = list.files("/ysm-gpfs/pi/zhao-data/zy92/GTEx_V8/GTEx_Analysis_v8_eQTL_covariates/")
idf = matrix(unlist(strsplit(cov_fl, "\\.")), ncol = 4, byrow = T)[,1]
##print(paste0("INFO: idf ", idf))
tissues = as.character(sapply(cov_fl, function(x) unlist(strsplit(x, '\\.'))[1]))

sink("summary_r2_all.sh")
for(t in tissues){
  cat("Rscript filter_gene_byfdr_all.R ",t,"\n")
}

sink()
