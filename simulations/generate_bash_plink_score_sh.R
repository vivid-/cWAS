sig <- read.table("/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/signature/LM22.txt",header=T,stringsAsFactors=F,sep="\t")
#sig <- readRDS("/ysm-gpfs/pi/zhao/wd262/sc/cWAS/sigmat/sigmat.rds")
genes <- colnames(sig)[-1]
#genes <- colnames(sig)
dat <- read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/gtexGene.txt",stringsAsFactors=F)
#c <- 0
tmp_df <- data.frame()
sink("remain_pred_expr_Whole_Blood_v8_1.sh")
#sink("pred_expr_Lung_ukbb.sh")
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
    if(file.exists(paste0(gene_fld,"Whole_Blood_weight.txt"))){
      if(!file.exists(gene_file)){
        tmp_df <- rbind(tmp_df,data.frame(chr=i,gene=g_1))
      }
    }
    if(!file.exists(gene_file)){
      next()
    }
    
    out_file <- paste0(gene_fld,"Whole_Blood_weight_ukbb_pred_gtex_v8")
    out_log <- paste0(gene_fld,"Whole_Blood_weight_ukbb_pred_gtex_v8.log") 
    if(file.exists(paste0(out_file,".profile"))){
       next()
    }
    #out_file <- paste0(gene_fld,"Whole_Blood_weight_ukbb_pred.expr.nopred")
 #   sent <- paste0("plink --bfile /ysm-gpfs/pi/zhao/yy496/ukb_imp_qc2/ukb_imp_qc2_chr",i," --score ",gene_file, " 1 4 2 \'header\' sum --out ",out_file,"; rm ",out_log)
    sent <- paste0("plink --bfile /ysm-gpfs/scratch60/wl382/GTEx/v8/genotype/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_maf0.01_grch38p7_chr",i," --score ",gene_file, " 1 4 2 sum --out ",out_file,"; rm ",out_log)
    #sent <- paste0("plink --bfile /ysm-gpfs/pi/zhao-data/wl382/cWAS/data/GTEX/v6/geno/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs --score ",gene_file, " 1 4 2 \'header\' sum --out ",out_file,"; rm ",out_log)
    cat(sent,"\n")
  }
}
sink()
write.table(tmp_df,file="remain_WB_genes_weights.txt",col.names=T,row.names=F,sep="\t",quote=F)
