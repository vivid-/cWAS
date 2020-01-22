args <- commandArgs()
signature_f <- args[6]
tissue <- args[7]

# load signature matrix
signature <- readRDS(file=signature_f)
genes <- colnames(signature)
# matrix for gene ids and gene symbols
gene_df<- read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/gtexGene.txt",stringsAsFactors=F)
# r2 summary
r2_df <- read.table(paste0("/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/model/r2/",tissue,"_all_summary.txt"),header=T,stringsAsFactors=F,sep="\t")
r2_df <- r2_df[which(p.adjust(r2_df$pval,method="fdr")<0.05),]
# select the gene_id_symbol df 
gene_df <- gene_df[match(genes,gene_df[,4]),]
gene_df_id <- gsub("\\.[0-9]+","",gene_df[,7])
r2_df_id <- gsub("\\.[0-9]+","",r2_df$gene)
# only focused on those signature genes with powerful imputation models
result_df <- gene_df[gene_df_id%in%r2_df_id,]

outputf <- paste0("pred_expr_v8_",tissue,".sh")
sink(outputf)
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
    gene_file <- paste0(gene_fld,tissue,"_weight_reformated.txt")
    #if(file.exists(paste0(gene_fld,"Whole_Blood_weight.txt"))){
    #  next()
    #}
    if(!file.exists(gene_file)){
      next()
    }
    
    out_file <- paste0(gene_fld,tissue,"_pred_gtex_v8")
    out_log <- paste0(gene_fld,tissue,"_pred_gtex_v8.log") 
    if(file.exists(paste0(out_file,".profile"))){
       next()
    }
    sent <- paste0("plink --bfile /ysm-gpfs/scratch60/wl382/GTEx/v8/genotype/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_maf0.01_grch38p7_chr",i," --score ",gene_file, " 1 4 2 sum --out ",out_file,"; rm ",out_log)
    cat(sent,"\n")
    }
}
sink()
