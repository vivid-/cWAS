#!/usr/bin/R

args <- commandArgs()
inputf <- args[6]
gene_df <- read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/gtexGene.txt",header=F,stringsAsFactors=F)
#colname_df <- read.table("/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/GTEX/v6/pred_expr/full_imputation_WB_txt_000",header=T,stringsAsFactors=F)
#colname_df <- read.table("summary_WB_ukbb_txt_000",header=T,stringsAsFactors=F)
colname_df <- read.table("summary_Lung_pred_ukbb_000",header=T,stringsAsFactors=F)
# read the imputed expression dataframe
#if(inputf != "summary_WB_ukbb_txt_000"){
if(inputf !="summary_Lung_pred_ukbb_000"){
  input <- read.table(inputf,header=F,stringsAsFactors=F)
  colnames(input) <- colnames(colname_df)
}else{
  input <- colname_df
}
#input <- read.table("summary_WB_ukbb.txt",header=T,stringsAsFactors=F)

# convert the colnames (ensembl id) to gene names using gene_df
new_colnames <- colnames(input)
new_colnames <- gsub("\\.[0-9]+","",new_colnames)
gene_df_names <- gsub("\\.[0-9]+","",gene_df[,7])
#gene_names <- gene_df[match(new_colnames[-1],gene_df[,7]),4] 
gene_names <- gene_df[match(new_colnames[-1],gene_df_names),4] 
new_colnames[-1] <- gene_names
colnames(input) <- new_colnames


# read the signature matrix
#signature <- read.table("/ysm-gpfs/pi/zhao-data/mc2792/sc_gwas/deconv/cybersort/LM22.txt",header=T,stringsAsFactors=F,sep="\t")
signature <- readRDS("/ysm-gpfs/pi/zhao/wd262/sc/cWAS/sigmat/sigmat.rds")
signature <- t(signature)
# extract the input matrix by gene names
#matched_ind <- match(signature$Gene.symbol,new_colnames)
matched_ind <- match(rownames(signature),new_colnames)
input_extract  <- input[,c(1,na.omit(matched_ind))]
# extract corresponding genes in signature matrix
#sig_extrac <- signature[match(colnames(input_extract)[-1],signature$Gene.symbol),]
sig_extrac <- signature[match(colnames(input_extract)[-1],rownames(signature)),]
cell_fracs <- data.frame()
# for each individual
for(i in 1:dim(input_extract)[1]){
  test <- data.frame(cbind(t(input_extract[i,-1]),sig_extrac[,-1]))
  colnames(test)[1] <- "y"
  a <- lm(y~.,data=test)
  cell_fracs <- rbind(cell_fracs,coef(a)[-1])
  colnames(cell_fracs) <- names(coef(a)[-1])
}

#cell_fracs[cell_fracs<0] <- 0
#a <- rowSums(cell_fracs)
#cell_fracs <- cell_fracs/a
cell_fracs$IID <- input$IID
outputf <- args[7]
write.table(cell_fracs,file=outputf,col.names=T,row.names=F,sep="\t",quote=F)
