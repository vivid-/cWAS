sig <- read.table("/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/signature/LM22.txt",header=T,stringsAsFactors=F,sep="\t")
genes <- colnames(sig)[-1]
dat <- read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/gtexGene.txt",stringsAsFactors=F)

#sink("pred_expr_lung_v8.sh")
result_df <- dat[match(genes,dat[,4]),]
total <- data.frame()
for(i in 1:22){
  tmp <- result_df[result_df[,1]==paste0("chr",i),]
  gs <- tmp[,7]
  gs <- na.omit(gs)
  #weight_fld <- paste0("/ysm-gpfs/scratch60/wl382/GTEx/v8/models/chr",i,"/")
  expr_fld <- paste0("/gpfs/loomis/scratch60/fas/zhao/zy92/GTEX/adjusted_expr1/chr",i,"/")
  fld_genes <- list.files(expr_fld)
  for(g in gs){
    g_1 <- gsub("\\.[0-9]+","",g)
    g_1 <- fld_genes[grep(g_1,fld_genes)]
    #cat(g_1,"\n") 
    gene_fld <- paste0(expr_fld,g_1,"/")
    gene_file <- paste0(gene_fld,"Whole_Blood.adj_expr")
    if(!file.exists(gene_file)){
      next()
    }
    tmp <- read.table(gene_file,header=F,stringsAsFactors=F)
    strsplit(tmp[,1],"-") -> n
    sapply(n,"[[",1) -> a
    sapply(n,"[[",2) -> b
    tmp[,1] <- paste0(a,"-",b) 
    if(length(total)==0){
      tmp -> total
      colnames(total) <- c("ID",result_df[which(result_df[,7]==g),4])
    }else{
      inter_ind <- intersect(tmp[,1],total[,1])
      total_o <- total[match(inter_ind,total[,1]),]
      cln_o <- colnames(total_o)
      tmp_1 <- tmp[match(inter_ind,tmp[,1]),2]
      total_o$new <- tmp_1
      colnames(total_o) <- c(cln_o,result_df[which(result_df[,7]==g),4])
      total <- total_o
    }
    #out_file <- paste0(gene_fld,"Whole_Blood_weight_ukbb_pred.expr")
    #out_log <- paste0(gene_fld,"Whole_Blood_weight_ukbb_pred.log") 
    #out_file <- paste0(gene_fld,"Whole_Blood_weight_ukbb_pred.expr.nopred")
    #sent <- paste0("plink --bfile /ysm-gpfs/pi/zhao/yy496/ukb_imp_qc2/ukb_imp_qc2_chr",i," --score ",gene_file, " 1 4 2 \'header\' sum --out ",out_file,"; rm ",out_log,";","rm ",out_file)
    #cat(sent,"\n")
  }
}
write.table(total,file="WB_v8_adj_exprs.txt",col.names=T,row.names=F,sep="\t",quote=F)
#sink()
