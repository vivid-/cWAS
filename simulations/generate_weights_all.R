sig <- read.table("/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/signature/LM22.txt",header=T,stringsAsFactors=F,sep="\t")
genes <- colnames(sig)[-1]
dat <- read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/gtexGene.txt",stringsAsFactors=F)
sel <- read.table("/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/model/r2/Whole_Blood_all_summary.txt",header=T,stringsAsFactors=F,sep="\t")
sel <- sel[which(p.adjust(sel$pval,method="fdr")<0.05),]
#sink("pred_expr_lung_v8.sh")
result_df <- dat[match(genes,dat[,4]),]
result_df_id <- gsub("\\.[0-9]+","",dat[,7])
sel_id <- gsub("\\.[0-9]","",sel$gene)
result_df <- result_df[result_df_id%in%sel_id,]
#total <- data.frame()
#genes_sel <- read.table("/ysm-gpfs/pi/zhao-data/wl382/cWAS/scripts/git_up/simulations/genes_mixed_sel_WB.txt",header=T,stringsAsFactors=F,sep="\t")
#dat_v6 <- read.table("/ysm-gpfs/project/wl382/snpPred_epi/Annotation/data/createDB/output/tables/Whole_Blood.elnt_allTissues_weight.txt",header=T,stringsAsFactors=F)
for(i in 1:22){
  total <- data.frame()
  tmp <- result_df[result_df[,1]==paste0("chr",i),]
  gs <- tmp[,7]
  gs <- na.omit(gs)
  #gs <- na.omit(genes_sel[genes_sel$chr==i,]$gene)
  weight_fld <- paste0("/ysm-gpfs/scratch60/wl382/GTEx/v8/models/chr",i,"/")
  fld_genes <- list.files(weight_fld)
  if(length(gs)==0){
    next()
  }
  for(g in gs){
    g_1 <- gsub("\\.[0-9]+","",g)
    g_1 <- fld_genes[grep(g_1,fld_genes)]
    #cat(g_1,"\n") 
    #if(genes_sel[genes_sel$gene==g,]$source=="v8"){
      gene_fld <- paste0(weight_fld,g_1,"/")
      gene_file <- paste0(gene_fld,"Whole_Blood_weight_reformated.txt")
      if(!file.exists(gene_file)){
        next()
      }
      tmp <- read.table(gene_file,header=F,stringsAsFactors=F,sep="\t")
      gene_symbol <- result_df[grep(g_1,result_df[,7]),4]
      cat(gene_symbol,"\n")
      colnames(tmp) <- c("rsID",gene_symbol,"ref","alt")
      tmp_1 <- data.frame(rsID=tmp$rsID,ref=tmp$ref,alt=tmp$alt,gene_symbol=tmp[,2])
    #}else{
    #  tmp <- unique(dat_v6[grep(g_1,dat_v6$gene),])

    #  gene_symbol <- result_df[grep(g_1,result_df[,7]),4]
    #  cat(gene_symbol,"\n")
    #  tmp_1 <- data.frame(rsID=tmp$rsid,ref=tmp$ref_allele,alt=tmp$eff_allele,gene_symbol=tmp$weight)
    #}
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

    #out_file <- paste0(gene_fld,"Whole_Blood_weight_ukbb_pred.expr")
    #out_log <- paste0(gene_fld,"Whole_Blood_weight_ukbb_pred.log") 
    #out_file <- paste0(gene_fld,"Whole_Blood_weight_ukbb_pred.expr.nopred")
    #sent <- paste0("plink --bfile /ysm-gpfs/pi/zhao/yy496/ukb_imp_qc2/ukb_imp_qc2_chr",i," --score ",gene_file, " 1 4 2 \'header\' sum --out ",out_file,"; rm ",out_log,";","rm ",out_file)
    #cat(sent,"\n")
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
  write.table(final,file=paste0("LM22_gtex_v8_weights_chr",i,"_sel.txt"),col.names=T,row.names=F,sep="\t",quote=F)
}
#final <- data.frame()
## match 
#identifers <- (unique(paste0(total[,1],"_",total[,2],"_",total[,3])))
#all <- paste0(total[,1],"_",total[,2],"_",total[,3])
#for(idens in identifers){
#   tmp <- total[which(identifers==idens),]
#   tmp_1 <- tmp[1,]
#   tmp_1[,-c(1:3)] <- colSums(tmp[,-c(1:3)])
#   final <- rbind(final,tmp_1)
#}

#write.table(final,"LM22_gtexv8_weights.txt",col.names=T,row.names=F,sep="\t",quote=F)
