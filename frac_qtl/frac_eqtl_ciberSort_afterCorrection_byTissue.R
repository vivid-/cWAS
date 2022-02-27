# calculate the association between cis SNPs and cell type fraction
library(data.table)
library("genio")
args <- commandArgs()
#chr <- args[6]
#gene <- args[7]
tissue <- args[6]

#gene_i <- gsub("\\.[0-9]+","",gene)
#cis_snp_dir <- paste0("/ysm-gpfs/project/wl382/GTEx_v8/genotype/cis_loc/",chr,"/",gene,"/",gene)
tissue_1 <- gsub("HCL_","",tissue)
tissue_1 <- gsub("_adult","",tissue_1)
tissue_1 <- gsub("_fetal","",tissue_1)
cis_snp_dir <- paste0("/gpfs/ysm/pi/zhao-data/wl382/cWAS/data/trans_qtl/",tissue_1)
frac_dir <- paste0("/gpfs/ysm/pi/zhao-data/wl382/cWAS/data/assayed_fractions/",tissue,"/")
frac_f <- paste0(frac_dir,tissue,"_adj_CiberSort_fracs.txt")

frac <- fread(frac_f,stringsAsFactors=F)
frac <- data.frame(frac)
geno <- read_plink(cis_snp_dir)
bed <- geno$X # snp by sample
bim <- data.frame(geno$bim)
fam <- data.frame(geno$fam)
# deal with the sample id in expr to match them with that in geno
a <- strsplit(frac[,1],"_")
sample_expr <- sapply(a,function(X) paste0(X[1],"-",X[2]))
inter_fam <- intersect(sample_expr,fam[,1])

# replace na with rowMeans
row_means <-  apply(bed,1,function(X) X[is.na(X)] <- mean(na.omit(X)))
inds <- which(is.na(bed),arr.ind=T)
bed[inds] <- row_means[inds[,1]]

#get the data for intersected samples
frac_1 <- frac[match(inter_fam,sample_expr),]
bed_1 <- bed[,match(inter_fam,fam[,1])]
fam_1 <- fam[match(inter_fam,fam[,1]),]
cell_num <- ncol(frac)-4
snp_num <- nrow(bed_1)


assoc_df <- data.frame()
for(i in 1:cell_num){
	for(j in 1:snp_num){
		tmp_md <- lm(frac_1[,i+1]~bed_1[j,])
		p <- summary(tmp_md)$coefficients[2,4]
		beta <- summary(tmp_md)$coefficients[2,1]
		se <- summary(tmp_md)$coefficients[2,2]
		tmp_df <- data.frame(snp = bim[j,],beta=beta,se=se,p=p,cell=colnames(frac)[i+1])
		assoc_df <- rbind(assoc_df,tmp_df)
	}
	#X_new = data.frame(y=frac[,i+1],X)
	#adj_frac = resid(lm(y~.,data=X_new))
	#frac_resid[,i+1] <- adj_frac
}

assoc_file <- paste0("/gpfs/ysm/pi/zhao-data/wl382/cWAS/data/frac_eqtl_HCL/",tissue,"cibersort_prop.eqtl")
write.table(assoc_df,file=assoc_file,col.names=T,row.names=F,sep="\t",quote=F)
