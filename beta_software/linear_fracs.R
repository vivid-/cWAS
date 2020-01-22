#!/usr/bin/R
args <- commandArgs()
signature_f <- args[6]
tissue <- args[7]

signature <- readRDS(file=signature_f)
genes <- colnames(signature)
# matrix for gene ids and gene symbols
gene_df <- read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/gtexGene.txt",stringsAsFactors=F)

# summarized pred_expr data
pred_expr <- read.table(paste0(tissue,"_summary_pred_gtex_v8.txt"),header=T,stringsAsFactors=F,sep="\t")
gene_ids <- gsub("\\.[0-9]+","",colnames(pred_expr)[-1])
gene_df_id <- gsub("\\.[0-9]+","",gene_df[,7])
gene_symbols <- gene_df[match(gene_ids,gene_df_id),4]
colnames(pred_expr)[-1] <- gene_symbols

# match the gene symbols with signature genes
sig_genes <- colnames(signature)
sig_cells <- rownames(signature)
signature <- t(signature)
matched_ind <- match(sig_genes,gene_symbols)
input_extract  <- pred_expr[,c(1,na.omit(matched_ind))]
sig_extrac <- signature[match(colnames(input_extract)[-1],sig_genes),]
cell_fracs <- data.frame()
# for each individual
for(i in 1:dim(input_extract)[1]){
  test <- cbind(t(input_extract[i,-1]),sig_extrac[,-1])
  colnames(test)[1] <- "y"
  a <- lm(y~.,data=test)
  cell_fracs <- rbind(cell_fracs,coef(a)[-1])
  colnames(cell_fracs) <- names(coef(a)[-1])
}

cell_fracs$FID <- pred_expr$IID
write.table(cell_fracs,file=paste0(tissue,"_frac_pred_expr_v8.txt"),col.names=T,row.names=F,sep="\t",quote=F)
