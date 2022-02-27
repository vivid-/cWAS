# cWAS: A Statistical Framework to Identify Cell Types Whose Genetically Regulated Proportions are Associated with Complex Diseases 




## cWAS
cWAS software (R library)
### Install required libraries
We are using `R version 3.5.0`  here. 
For R, `genio` and `MASS` library are required, which can be installed via 
``` R
> install.packages("genio")
> install.packages("R.utils")
```
### Install the cWAS library from source
```R
> install.packages("cWAS.tar.gz",repos=NULL,type="source")
```
### Formats of input data
1. gene expression imputation weights data (weight_f)···

   The file should have five columns of data included (chr id	posg	pos	ref	alt	weight	gene	tissue) as shown below:
   ```R
   chr	id	posg	pos	ref	alt	weight	gene	tissue
   chr1	rs12071538	0	22674247	T	C	0.025650281962564	ENSG00000007968.6	Whole_Blood
   chr1	rs115853779	0	156236904	G	A	-0.0590590453038043	ENSG00000027644.4	Whole_Blood
   chr1	rs115617321	0	156244942	G	A	0.0228122359672935	ENSG00000027644.4	Whole_Blood
   chr1	rs201707266	0	156800790	TG	T	0.0493494890914809	ENSG00000027644.4	Whole_Blood
   ```
   
2. signature gene expression matrix (signature_f)

   RDS file includes a dataframe with column names being gene symbols and row names being cell types.
   
3. the file with gene symbols and corresponding ensembl ids (gene_df_f)
   ```R
   #chrom	chromStart	chromEnd	name	score	strand	geneId	geneType	expCount	expScores
   chr1	11868	13052	DDX11L1	15	+	ENSG00000223972.4	transcribed_unprocessed_pseudogene	53	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.239,0,0,0,0,
   chr1	14969	29806	WASH7P	427	-	ENSG00000227232.4	unprocessed_pseudogene	53	   6.886,6.083,4.729,5.91,6.371,6.007,8.768,4.202,4.455,4.64,10.097,10.619,6.108,5.037,5.018,4.808,4.543,4.495,5.576,4.57,8.275,4.707,2.55,9.091,9.885,8.17,7.392,7.735,5.353,7.124,8.617,3.426,2.375,7.669,3.826,7.094,6.365,3.263,10.723,10.507,4.843,9.193,13.25,11.635,11.771,8.641,10.448,6.522,9.313,10.304,9.987,9.067,6.12,
   ```
   
 4. gene expression imputation peformance file (r2_f)
    ```R
    gene	tissue	rsq	pval	chr
    ENSG00000000460.16	Whole_Blood	0.126802007448886	1.84129251898404e-21	1
    ENSG00000001460.17	Whole_Blood	0.122792740236351	8.63833212300849e-21	1
    ```
 5. gwas weight file (gwas_f)
 
    It should include the following five columns, where A2 is the non-effect allele
    ```R
    SNP A1 A2 Z N P
    rs10875231 T G -0.156 63926 0.8760330212212881
    rs6678176 T C 0.053 63926 0.9577319077150701
    rs78286437 T C -0.414 63926 0.6788741068331019
    rs144406489 A G -0.1 63926 0.920344325445942
    ```
  6. reference genome location/bim file (ref_geno_f)
  
     It should have the same data format as that of plink bim file 


## simulations
Code for generating simulated phenotypes and also GWAS summary stats

## build_model
Code for generating gene expression imputation weights based on GTEx v8 data
