sig <- read.table("/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/signature/LM22.txt",header=T,stringsAsFactors=F,sep="\t")
genes <- colnames(sig)[-1]
dat <- read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/gtexGene.txt",stringsAsFactors=F)

sink("pred_expr_lung_v8_cp.sh")
result_df <- dat[match(genes,dat[,4]),]
for(i in 1:22){
  tmp <- result_df[result_df[,1]==paste0("chr",i),]
  gs <- tmp[,7]
  gs <- na.omit(gs)
  weight_fld <- paste0("/ysm-gpfs/scratch60/wl382/GTEx/v8/models/chr",i,"/")
  fld_genes <- list.files(weight_fld)
  for(g in gs){
    g_1 <- gsub("\\.[0-9]+","",g)
    g_1 <- fld_genes[grep(g_1,fld_genes)]
    #cat(g_1,"\n") 
    gene_fld <- paste0(weight_fld,g_1,"/")
    gene_file <- paste0(gene_fld,"Whole_Blood_weight_reformated.txt")
    result_file <- paste0(gene_fld,"Whole_Blood_weight_reformated_simul.txt")
    if(!file.exists(gene_file)){
      next()
    }
    out_file <- paste0(gene_fld,"Whole_Blood_weight_ukbb_pred.expr")
    out_log <- paste0(gene_fld,"Whole_Blood_weight_ukbb_pred.log") 
    out_file <- paste0(gene_fld,"Whole_Blood_weight_ukbb_pred.expr.nopred")
#    sent <- paste0("plink --bfile /ysm-gpfs/pi/zhao/yy496/ukb_imp_qc2/ukb_imp_qc2_chr",i," --score ",gene_file, " 1 4 2 \'header\' sum --out ",out_file,"; rm ",out_log,";","rm ",out_file)
    sent <- paste0("cp ",gene_file," ",result_file)
    cat(sent,"\n")
  }
}
sink()
