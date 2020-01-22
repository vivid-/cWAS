#!/usr/bin/R
args <- commandArgs()
signature_f <- args[6]
tissue <- args[7]

# load signature matrix
signature <- readRDS(file=signature_f)
genes <- colnames(signature)
# matrix for gene ids and gene symbols
gene_df <- read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/gtexGene.txt",stringsAsFactors=F)
# r2 summary
r2_df <- read.table(paste0("/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/model/r2/",tissue,"_all_summary.txt"),header=T,stringsAsFactors=F,sep="\t")
r2_df <- r2_df[which(p.adjust(r2_df$pval,method="fdr")<0.05),]
# select the gene_id_symbol df 
gene_df <- gene_df[match(genes,gene_df[,4]),]
gene_df_id <- gsub("\\.[0-9]+","",gene_df[,7])
r2_df_id <- gsub("\\.[0-9]+","",r2_df$gene)
# only focused on those signature genes with powerful imputation models
result_df <- gene_df[gene_df_id%in%r2_df_id,]


all_genes <- data.frame()
for(i in 1:22){
  tmp <- result_df[result_df[,1]==paste0("chr",i),]
  gs <- tmp[,7]
  gs <- na.omit(gs)
  weight_fld <- paste0("/ysm-gpfs/scratch60/wl382/GTEx/v8/models/chr",i,"/")
  fld_genes <- list.files(weight_fld)
  for(g in gs){
    g_1 <- gsub("\\.[0-9]+","",g)
    if(length(grep(g_1,sel$gene))==0){
       next()
    }
    g_1 <- fld_genes[grep(g_1,fld_genes)]
    all_genes <- rbind(all_genes,data.frame(chr=i,gene=g_1))

    out_file <- paste0(gene_fld,tissue,"_pred_gtex_v8.profile")
    if(!file.exists(out_file)){
       next()
    }
    tmp<- read.table(out_file,header=T,stringsAsFactors=F)
    if(length(result)==0){
      result <- tmp[,c("IID","SCORESUM")]
      colnames(result) <- c("IID",g_1)
    }else{
      names <- colnames(result)
      result <- cbind(result,tmp[,"SCORESUM"])
      colnames(result) <- c(names,g_1)
    }
  }
}

write.table(result,file=paste0(tissue,"_summary_pred_gtex_v8.txt"),col.names=T,row.names=F,sep="\t",quote=F)

