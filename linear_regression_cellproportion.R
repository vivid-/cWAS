#!/usr/bin/R

args <- commandArgs()
inputf <- args[6]
gene_df <- read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/gtexGene.txt",header=F,stringsAsFactors=F)
colname_df <- read.table("/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/GTEX/v6/pred_expr/full_imputation_WB_txt_000",header=T,stringsAsFactors=F)

# read the imputed expression dataframe
if(inputf != "/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/GTEX/v6/pred_expr/full_imputation_WB_txt_000"){
  input <- read.table(inputf,header=F,stringsAsFactors=F)[,-1]
  colnames(input) <- colnames(colname_df)
}else{
  input <- colname_df
}

# convert the colnames (ensembl id) to gene names using gene_df
new_colnames <- colnames(input)
gene_names <- gene_df[match(new_colnames[-1],gene_df[,7]),4] 
new_colnames[-1] <- gene_names
colnames(input) <- new_colnames


# read the signature matrix
signature <- read.table("/ysm-gpfs/pi/zhao/mc2792/sc_gwas/deconv/cybersort/LM22.txt",header=T,stringsAsFactors=F,sep="\t")
# extract the input matrix by gene names
matched_ind <- match(signature$Gene.symbol,new_colnames)
input_extract  <- input[,c(1,na.omit(matched_ind))]
# extract corresponding genes in signature matrix
sig_extrac <- signature[match(colnames(input_extract)[-1],signature$Gene.symbol),]

cell_fracs <- data.frame()
# for each individual
for(i in 1:dim(input_extract)[1]){
  test <- cbind(t(input_extract[i,-1]),sig_extrac[,-1])
  colnames(test)[1] <- "y"
  a <- lm(y~.,data=test)
  cell_fracs <- rbind(cell_fracs,coef(a)[-1])
  colnames(cell_fracs) <- names(coef(a)[-1])
}

#cell_fracs[cell_fracs<0] <- 0
#a <- rowSums(cell_fracs)
#cell_fracs <- cell_fracs/a
cell_fracs$FID <- input$FID
outputf <- args[7]
write.table(cell_fracs,file=outputf,col.names=T,row.names=F,sep="\t",quote=F)
