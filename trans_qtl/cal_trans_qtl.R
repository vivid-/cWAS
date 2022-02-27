#!/usr/bin/R
# calculate the eQTL effect size and p values for each gene-SNP pair
library(data.table)
library("genio")
args <- commandArgs()
chr <- args[6]
gene <- args[7]
tissue <- args[8]
gene_i <- gsub("\\.[0-9]+","",gene)
cis_snp_dir <- "/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/trans_qtl/"
cis_snp_f <- paste0(cis_snp_dir,tissue)

expr_dir <- paste0("/gpfs/loomis/scratch60/zhao/wl382/GTEx_V8/adjusted_expr/",chr,"/")

gene_i <- list.files(expr_dir,pattern=gene_i)

expr_f <- paste0(expr_dir,"/",gene_i,"/",tissue,".adj_expr")
if(file.exists(expr_f)){
  expr_dir <- paste0("/gpfs/loomis/scratch60/zhao/ky292/GTEx_V8/adjusted_expr/",chr,"/")

  gene_i <- list.files(expr_dir,pattern=gene_i)

  expr_f <- paste0(expr_dir,"/",gene_i,"/",tissue,".adj_expr")
}



geno <- read_plink(cis_snp_f)
expr <- read.table(expr_f,header=F,stringsAsFactors=F)

bed <- geno$X # snp by sample
bim <- data.frame(geno$bim)
fam <- data.frame(geno$fam)

# deal with the sample id in expr to match them with that in geno
a <- strsplit(expr[,1],"-")
sample_expr <- sapply(a,function(X) paste0(X[1],"-",X[2]))
inter_fam <- intersect(sample_expr,fam[,1])

# replace na with rowMeans
row_means <-  apply(bed,1,function(X) X[is.na(X)] <- mean(na.omit(X)))
inds <- which(is.na(bed),arr.ind=T)
bed[inds] <- row_means[inds[,1]]

#get the data for intersected samples
expr_1 <- expr[match(inter_fam,sample_expr),]
bed_1 <- bed[,match(inter_fam,fam[,1])]
fam_1 <- fam[match(inter_fam,fam[,1]),]

all_df <- data.frame()
for(i in 1:nrow(bed_1)){
	tmp_md <- lm(expr_1[,2]~bed_1[i,])
        
	p <- summary(tmp_md)$coefficients[2,4]
    beta <- summary(tmp_md)$coefficients[2,1]
    se <- summary(tmp_md)$coefficients[2,2]
    tmp_df <- data.frame(snp = bim[i,2],beta=beta,se=se,p=p,gene=gene)
    all_df <- rbind(all_df,tmp_df)
}

# get the permuted p values
y = expr_1[,2]
reps = 1000000
k = length(y)
expr_1_permute <- unique(t(sapply(1:reps, function(x) sample(y, k))))
permutes_p <- c()


out_dir <- paste0("/ysm-gpfs/scratch60/ky292/cWAS/data/trans_qtl/",chr,"/",gene)
dir.create(out_dir ,showWarnings = FALSE,recursive=TRUE)
outf <- paste0("/ysm-gpfs/scratch60/ky292/cWAS/data/trans_qtl/",chr,"/",gene,"/",tissue,"_trans_eqtl.txt")

write.table(all_df,file=outf,col.names=T,row.names=F,sep="\t",quote=F)


