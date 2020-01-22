#!/usr/bin/R

args <- commandArgs()
chr <- args[6]
sig <- read.table("/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/signature/LM22.txt",header=T,stringsAsFactors=F,sep="\t")
genes_sel <- colnames(sig)[-1]
genes <- read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/gtexGene.txt",stringsAsFactors=F)

dat <- read.table("/ysm-gpfs/project/wl382/snpPred_epi/Annotation/data/createDB/output/tables/Whole_Blood.elnt_allTissues_weight.txt",header=T,stringsAsFactors=F)

#genes <- read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/gtexGene.txt",header=F,stringsAsFactors=F,sep="\t")

sel <- read.table("/ysm-gpfs/pi/zhao/wl382/snpPred_epi/Annotation/data/model_r2/selected_models_elnt_Whole_Blood.txt",header=T,stringsAsFactors=F,sep="\t")
dat <- dat[dat$gene %in%sel$gene,]
#genes <- genes[genes[,1]==chr,]
#gene_names <- genes[,7]
all_genes <- data.frame()
for(chr in 1:22){
  genes_1 <- genes[genes[,1]==paste0("chr",chr),]
  result_df <- genes_1[na.omit(match(genes_sel,genes_1[,4])),]
  if(length(na.omit(match(genes_sel,genes_1[,4])))==0){
    next()
  }
  total <- data.frame()
  for(i in 1:nrow(result_df)){
    g <- result_df[i,7]
    g_0 <- gsub("\\.[0-9]+","",g)
    if(length(grep(g_0,dat$gene))==0){
      next()
    }
    all_genes <- rbind(all_genes,data.frame(chr=chr,gene=g))
   # total <- data.frame()
     tmp <- dat[grep(g_0,dat$gene),]
     gene_symbol <- result_df[grep(g,result_df[,7]),4]
     cat(gene_symbol,"\n")
     #colnames(tmp) <- c("rsID",gene_symbol,"ref","alt")
     tmp_1 <- data.frame(rsID=tmp$rsid,ref=tmp$ref_allele,alt=tmp$eff_allele,gene_symbol=tmp$weight)
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
  write.table(final,file=paste0("LM22_gtexv6_weights_sel_chr",chr,".txt"),col.names=T,row.names=F,sep="\t",quote=F)

}
write.table(all_genes,file="WB_v6_genes_sel.txt",col.names=T,row.names=F,sep="\t",quote=F)
#for(g in gene_names){
#  if(!g %in% dat$gene){
#     next()
#  }
#  if(!g %in% sel$gene){
#     next()
#  }
#  tmp <- unique(dat[dat$gene==g,])
#  file_out <- paste0("/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/GTEX/v6/weights/Whole_Blood/",chr,"/",g,"_weights.txt")
#  write.table(tmp,file=file_out,col.names=F,row.names=F,sep="\t",quote=F)
#}

