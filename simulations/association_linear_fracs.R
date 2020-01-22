#!/usr/bin/R

# read frac files first
#files <- list.files("/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/GTEX/v6/fracs_linearRegression/",pattern="full_imputation_WB_cell_fracs_")
#files <- list.files("/ysm-gpfs/pi/zhao-data/wl382/cWAS/scripts/git_up/simulations/",pattern="frac_WB_ukbb_txt_")
files <- list.files("/ysm-gpfs/pi/zhao-data/wl382/cWAS/scripts/git_up/simulations/",pattern="frac_Lung_pred_ukbb_")
result_df <- data.frame()
for(file in files){
  tmp_file <- paste0("/ysm-gpfs/pi/zhao-data/wl382/cWAS/scripts/git_up/simulations/",file)
  tmp <- read.table(tmp_file,header=T,stringsAsFactors=F)

  result_df <- rbind(result_df,tmp)
}
args <- commandArgs()
h2 <- as.numeric(args[6])
# sample size  
n <- as.numeric(args[7])
outputf <- args[8]
#sel_sample <- sample(1:nrow(result_df),n)
fam <- read.table("/ysm-gpfs/project/wl382/cWAS/data/simul/geno/h2_0.1_all.fam",header=F,stringsAsFactors=F)
sel_sample <- match(fam[,1],result_df$IID)
a <- !is.na(result_df)
# select those cell types/columns having most individuals having non-NA fractions
sum_ind <- colSums(a)
cell_inds <- which(sum_ind==max(sum_ind))
result_df_1 <- result_df[,cell_inds]
# select the cell types with the largest variance
#a <- apply(result_df_1[,-ncol(result_df_1)],2,var)
#cell_ind <- which(a==max(a))
#cell_ind <- which(colnames(result_df_1)=="Macrophages.M1")
cell_ind <- which(colnames(result_df_1)=="Basal")
# treat the cell fracs as predictors
cat(cell_ind,"\n")
x <- result_df_1[sel_sample,cell_ind]

# sample beta for cell fraction and repeat multiple times
for(rep in 1:300){
  if(h2!=0){
  beta <- rnorm(1)
  xb <- x*beta
  sigma_e <- (1-h2)/h2*var(xb)
  e <- rnorm(n,mean=0,sd=sqrt(sigma_e))
  pheno <- xb+e
  }
  else{
    pheno <- rnorm(1)
    beta = 0
  }
  result_df <- data.frame(FID=result_df_1[sel_sample,]$IID,IID=result_df_1[sel_sample,]$IID,pheno=pheno)
  write.table(result_df,file=paste0(outputf,"_",rep,".pheno"),col.names=T,row.names=F,sep="\t",quote=F)   
  param <- data.frame(cell=colnames(result_df_1)[cell_ind],beta=beta)
  write.table(param,file=paste0(outputf,"_",rep,".param"),col.names=T,row.names=F,sep="\t",quote=F)
}

