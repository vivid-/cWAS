sig <- read.table("/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/signature/LM22.txt",header=T,stringsAsFactors=F,sep="\t")
genes <- colnames(sig)[-1]
#sig <- readRDS("/ysm-gpfs/pi/zhao/wd262/sc/cWAS/sigmat/sigmat.rds")
#genes <- colnames(sig)[-1]
#genes <- colnames(sig)
dat <- read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/gtexGene.txt",stringsAsFactors=F)
#sel <- read.table("/ysm-gpfs/pi/zhao/wl382/snpPred_epi/Annotation/data/model_r2/selected_models_elnt_Whole_Blood.txt",header=T,stringsAsFactors=F,sep="\t")


result <- data.frame()
#sink("pred_expr_lung_v8.sh")
result_df <- dat[match(genes,dat[,4]),]
sel <- read.table("/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/model/r2/Whole_Blood_all_summary.txt",header=T,stringsAsFactors=F,sep="\t")
sel <- sel[which(p.adjust(sel$pval,method="fdr")<0.05),]
#sink("pred_expr_lung_v8.sh")
#result_df <- dat[match(genes,dat[,4]),]
result_df_id <- gsub("\\.[0-9]+","",result_df[,7])
sel_id <- gsub("\\.[0-9]","",sel$gene)
result_df <- result_df[result_df_id%in%sel_id,]
gene_num <- 0
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
    gene_num <- gene_num+1
    g_1 <- fld_genes[grep(g_1,fld_genes)]
    all_genes <- rbind(all_genes,data.frame(chr=i,gene=g_1))
    #cat(g_1,"\n") 
    gene_fld <- paste0(weight_fld,g_1,"/")
    #gene_file <- paste0(gene_fld,"Whole_Blood_weight_reformated.txt")
    #if(!file.exists(gene_file)){
    #  next()
    #}
    out_file <- paste0(gene_fld,"Whole_Blood_weight_ukbb_pred_gtex_v8.profile")
    #out_file <- paste0("Whole_Blood_weight_ukbb_pred.expr.nopred.profile")
    #out_file <- paste0(gene_fld,"Whole_Blood_weight_ukbb_pred.expr.nopred.profile")
    #out_file <- paste0(gene_fld,"Lung_weight_ukbb_pred_gtex_v8.profile")
    #out_file <- paste0(gene_fld,"Lung_weight_ukbb_pred.profile")    
    if(!file.exists(out_file)){
       next()
    }
    #out_log <- paste0(gene_fld,"Whole_Blood_weight_ukbb_pred.log") 
    #out_file <- paste0(gene_fld,"Whole_Blood_weight_ukbb_pred.expr.nopred")
    #sent <- paste0("plink --bfile /ysm-gpfs/pi/zhao/yy496/ukb_imp_qc2/ukb_imp_qc2_chr",i," --score ",gene_file, " 1 4 2 \'header\' sum --out ",out_file,"; rm ",out_log,";","rm ",out_file)
    #cat(sent,"\n")
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
cat(gene_num,"\n")
#sink()
write.table(all_genes,file="genes_v8_sel_WB.txt",col.names=T,row.names=F,sep="\t",quote=F)
write.table(result,file="summary_WB_pred_v8_sel.txt",col.names=T,row.names=F,sep="\t",quote=F)
