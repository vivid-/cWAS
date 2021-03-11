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



## simulations
Code for generating simulated phenotypes and also GWAS summary stats

## build_model
Code for generating gene expression imputation weights based on GTEx v8 data
