#!/usr/bin/R

# read frac files first
files <- list.files("/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/GTEX/v6/fracs_linearRegression/",pattern="full_imputation_WB_cell_fracs_")
result_df <- data.frame()
for(file in files){
  tmp_file <- paste0("/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/GTEX/v6/fracs_linearRegression/",file)
  tmp <- read.table(tmp_file,header=T,stringsAsFactors=F)

  result_df <- rbind(result_df,tmp)
}
args <- commandArgs()
# read input phenotype files
pheno_file <- args[6]
pheno <- read.table(pheno_file,header=F,stringsAsFactors=F)

# match samples in phenotype file and cell fraction file
pheno_matched <- pheno[match(result_df$FID,pheno[,1]),6]
result_df$pheno <- pheno_matched

assoc_df <- data.frame()
for(i in 1:22){
  tmp <- result_df[,i]
  tmp_model <- lm(pheno_matched~tmp)
  gamma <- summary(tmp_model)$coefficients[2,1]
  p <- summary(tmp_model)$coefficients[2,4]
  tmp_df <- data.frame(cell=colnames(result_df)[i],gamma=gamma,p=p)
  assoc_df <- rbind(assoc_df,tmp_df)
}
assoc_df_sig <- assoc_df[assoc_df$p<0.05/dim(assoc_df)[1],]

#tmp <- lm(pheno~.,data=result_df)
outputf_file <- args[7]
#save(tmp,file=outputf_file)
write.table(assoc_df,file=outputf_file,col.names=T,row.names=F,sep="\t",quote=F)
