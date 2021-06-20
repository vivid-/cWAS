\name{cWAS}
\alias{cWAS}
\title{cWAS:A Statistical Framework to Identify Cell Types Whose Genetically Regulated Proportions are Associated with Complex Diseases}
\usage{
cWAS(gene_df_f,weight_f,r2_f,signature_f, gwas_f,ref_geno_f,output_f)
cWAS_test(ref_geno_f,S_f,Z_f,betaS_f,out_f)
}

\arguments{
\item{gene_df_f}{file location for the file including gene ensembl ids and gene ids}
\item{weight_f}{file location for imputation weights}
\item{r2_f}{file location for gene expression imputation R2}
\item{signature_rds_f}{file location for signature gene expression matrix rds file}
\item{gwas_f}{file location for GWAS summary stats}
\item{ref_geno_f}{file location for reference genome}
\item{output_f}{output file prefix after the SNP and gene matching process}
\item{S_f}{file location for signature matrix from the output of cWAS()}
\item{Z_f}{file location for GWAS summary stats from the output of cWAS()}
\item{betaS_f}{file location for betaS from the output of cWAS()}
\item{out_f}{file location for output association test results}
\item{permute}{logic variable indicating whether the permutation-like test should be conducted or not, the default input is FALSE}
\item{n_permute}{the number of permutation times, the default input is 1e6}
}
\description{
cWAS() is to match SNPs and genes in imputation weights, signature gene expression matrix, reference genotype input and GWAS summary stats
cWAS_test() is used to test the association between genetic-regulated cell type proportions and the phenotype of interests using the output from cWAS()
}
\examples{
cWAS(gene_df_f,weight_f,r2_f,signature_f, gwas_f,ref_geno_f,output_f)
S_f = paste0(output_f,".signature")
Z_f = paste0(output_f,".gwas")
betaS_f = paste0(output_f,".betaS")
# generating the permutation-like test results as well
cWAS_test(ref_geno_f,S_f,Z_f,betaS_f,out_f,permute=TRUE,n_permute=1e6)
}
