rep=$1
h2=$2
module purge
mkdir -p /ysm-gpfs/project/wl382/cWAS/data/simul/gwas/h2_${h2}/${rep}/
for i in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,22};do
Rscript /ysm-gpfs/pi/zhao-data/wl382/cWAS/scripts/git_up/beta_software/match_gwas_expr_ref.R --out-expr-beta=/ysm-gpfs/project/wl382/cWAS/data/simul/gwas/h2_${h2}/${rep}/expr_beta_v8_chr${i}.txt --out-gwas-beta=/ysm-gpfs/project/wl382/cWAS/data/simul/gwas/h2_${h2}/${rep}/gwas_chr${i}.txt --expr-beta=/ysm-gpfs/pi/zhao-data/wl382/cWAS/scripts/git_up/simulations/LM22_gtexv8_weights_chr${i}.txt --signature-matrix=/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/signature/LM22.txt --gwas-beta=/ysm-gpfs/project/wl382/cWAS/data/simul/gwas/h2_${h2}/${rep}.gwas.txt  --ref-geno=/ysm-gpfs/scratch60/wl382/GTEx/v8/genotype/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_maf0.01_grch38p7 --out-sig=/ysm-gpfs/project/wl382/cWAS/data/simul/gwas/h2_${h2}/${rep}/LM22_chr${i}.txt --out-ref-geno=/ysm-gpfs/project/wl382/cWAS/data/simul/gwas/h2_${h2}/${rep}/gtex_v8_matched_chr${i}.bim
done

# get expr_beta_S
for i in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,22};do
Rscript /ysm-gpfs/pi/zhao-data/wl382/cWAS/scripts/git_up/beta_software/get_expr_beta_s.R --out=/ysm-gpfs/project/wl382/cWAS/data/simul/gwas/h2_${h2}/${rep}/expr_beta_S_chr${i}.txt --expr-beta=/ysm-gpfs/project/wl382/cWAS/data/simul/gwas/h2_${h2}/${rep}/expr_beta_v8_chr${i}.txt --signature-matrix=/ysm-gpfs/project/wl382/cWAS/data/simul/gwas/h2_${h2}/${rep}/LM22_chr${i}.txt
done

# summary expr_beta 
Rscript /ysm-gpfs/pi/zhao-data/wl382/cWAS/scripts/git_up/beta_software/summary_expr_beta.R /ysm-gpfs/project/wl382/cWAS/data/simul/gwas/h2_${h2}/${rep}/
Rscript /ysm-gpfs/pi/zhao-data/wl382/cWAS/scripts/git_up/beta_software/summary_beta_S.R /ysm-gpfs/project/wl382/cWAS/data/simul/gwas/h2_${h2}/${rep}/
#script summary_files.R /cWAS/data/simul/gwas/${rep}/
Rscript /ysm-gpfs/pi/zhao-data/wl382/cWAS/scripts/git_up/beta_software/summary_files.R /ysm-gpfs/project/wl382/cWAS/data/simul/gwas/h2_${h2}/${rep}/
# plink
plink --bfile /ysm-gpfs/scratch60/wl382/GTEx/v8/genotype/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_maf0.01_grch38p7 --extract /ysm-gpfs/project/wl382/cWAS/data/simul/gwas/h2_${h2}/${rep}/gtex_v8_all.bim --make-bed --out /ysm-gpfs/project/wl382/cWAS/data/simul/gwas/h2_${h2}/${rep}/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs_all
module load R/3.6.1-foss-2018b;Rscript /ysm-gpfs/pi/zhao-data/wl382/cWAS/scripts/git_up/beta_software/assoc_test_v2.R --out=/ysm-gpfs/project/wl382/cWAS/data/simul/gwas/h2_${h2}/${rep}/assoc_LM22 --signature-matrix=/ysm-gpfs/project/wl382/cWAS/data/simul/gwas/h2_${h2}/${rep}/LM22_new_all.txt --gwas-beta=/ysm-gpfs/project/wl382/cWAS/data/simul/gwas/h2_${h2}/${rep}/gwas_all.txt --ref-geno=/ysm-gpfs/project/wl382/cWAS/data/simul/gwas/h2_${h2}/${rep}/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs_all --ref-cellFraction=/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/GTEX/v6/ref_cell_fracs_predExpr_Whole_Blood.txt --tol-inv=0.01 --betaS-file=/ysm-gpfs/project/wl382/cWAS/data/simul/gwas/h2_${h2}/${rep}/expr_beta_S_all.txt;rm /ysm-gpfs/project/wl382/cWAS/data/simul/gwas/h2_${h2}/${rep}/*bim;rm /ysm-gpfs/project/wl382/cWAS/data/simul/gwas/h2_${h2}/${rep}/*bed;rm /ysm-gpfs/project/wl382/cWAS/data/simul/gwas/h2_${h2}/${rep}/*fam;rm -rf /ysm-gpfs/project/wl382/cWAS/data/simul/gwas/h2_${h2}/${rep}/*chr*txt
