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

# for each chromosome
for(i in 1:22){
  total <- data.frame()
  tmp <- result_df[result_df[,1]==paste0("chr",i),]
  gs <- tmp[,7]
  gs <- na.omit(gs)
  weight_fld <- paste0("/ysm-gpfs/scratch60/wl382/GTEx/v8/models/chr",i,"/")
  fld_genes <- list.files(weight_fld)
  if(length(gs)==0){
    next()
  }
  
  for(g in gs){
    g_1 <- gsub("\\.[0-9]+","",g)
    g_1 <- fld_genes[grep(g_1,fld_genes)]
    gene_fld <- paste0(weight_fld,g_1,"/")
    gene_file <- paste0(gene_fld,"Whole_Blood_weight_reformated.txt")
    if(!file.exists(gene_file)){
      next()
    }
    tmp <- read.table(gene_file,header=F,stringsAsFactors=F,sep="\t")
    gene_symbol <- result_df[grep(g_1,result_df[,7]),4]
    colnames(tmp) <- c("rsID",gene_symbol,"ref","alt")
    tmp_1 <- data.frame(rsID=tmp$rsID,ref=tmp$ref,alt=tmp$alt,gene_symbol=tmp[,2])
    colnames(tmp_1)[4] <- gene_symbol
    if(length(total)==0){
      total <- tmp_1
    }else{
      snp_df <- tmp_1[,c("rsID","ref","alt")]
      gene_existing <- colnames(total)[!colnames(total) %in% c("rsID","ref","alt")]
      mm <- matrix(0,dim(snp_df)[1], length(gene_existing))
      snp_df <- cbind(snp_df,data.frame(mm))
      colnames(snp_df) <- c("rsID","ref","alt",gene_existing)
      gene_adding <- colnames(tmp_1)[!colnames(tmp_1) %in% c("rsID","ref","alt")]
      ur <- matrix(0,dim(total)[1], length(gene_adding))
      
      snp_df <- cbind(snp_df,tmp_1[,!colnames(tmp_1) %in% c("rsID","ref","alt")])
      colnames(snp_df) <- c("rsID","ref","alt",gene_existing,gene_adding)
      coln <- colnames(total)

      total <- cbind(total,data.frame(ur))
      colnames(total) <- c(coln,gene_adding)
      total <- rbind(total,snp_df)
    }
    
    if(length(total)==0){
      next()
    }
    final <- data.frame()
    identifers <- unique(paste0(total[,1],"_",total[,2],"_",total[,3]))
    all <- paste0(total[,1],"_",total[,2],"_",total[,3])
    for(idens in identifers){
      tmp <- total[which(all==idens),]
      tmp_1 <- tmp[1,]
      if(nrow(tmp)!=1){
        tmp_1[,-c(1:3)] <- colSums(tmp[,-c(1:3)])
      }
      final <- rbind(final,tmp_1)
    }
  
    write.table(final,file=paste0(tissue,"_gtex_v8_weights_chr",i,"_sel.txt"),col.names=T,row.names=F,sep="\t",quote=F)
  }
}
