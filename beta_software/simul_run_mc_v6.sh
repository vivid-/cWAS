rep=$1 #replication
h2=$2 #disease

# original /ysm-gpfs/project/wl382/cWAS/data/simul/gwas
path1=/home/mc2792/scratch60/wei/analysis 

module purge
module load R/3.6.1-foss-2018b

mkdir -p $path1/Lung_h2_${h2}/${rep}/

for i in {1..22}
do
    if test -f "LM22_chr$i.txt"
        then 
            echo "chr $i$ done!"
        else 
            Rscript /ysm-gpfs/pi/zhao-data/wl382/cWAS/scripts/git_up/beta_software/match_gwas_expr_ref.R --out-expr-beta=$path1/Lung_h2_${h2}/${rep}/expr_beta_v8_chr${i}.txt --out-gwas-beta=$path1/Lung_h2_${h2}/${rep}/gwas_chr${i}.txt --expr-beta=/ysm-gpfs/pi/zhao-data/wl382/cWAS/scripts/git_up/simulations/LM22_gtexv6_weights_chr${i}.txt --signature-matrix=/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/signature/LM22.txt --gwas-beta=$path1/Lung_h2_${h2}/${rep}.gwas.txt --ref-geno=/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/GTEX/v6/geno/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs_hg19 --out-sig=$path1/Lung_h2_${h2}/${rep}/LM22_chr${i}.txt --out-ref-geno=$path1/Lung_h2_${h2}/${rep}/gtex_v8_matched_chr${i}.bim
    fi
done

# get expr_beta_S
for i in {1..22};do
Rscript /ysm-gpfs/pi/zhao-data/wl382/cWAS/scripts/git_up/beta_software/get_expr_beta_s.R --out=$path1/Lung_h2_${h2}/${rep}/expr_beta_S_chr${i}.txt --expr-beta=$path1/Lung_h2_${h2}/${rep}/expr_beta_v8_chr${i}.txt --signature-matrix=$path1/Lung_h2_${h2}/${rep}/LM22_chr${i}.txt
done

# summary expr_beta 
Rscript /ysm-gpfs/pi/zhao-data/wl382/cWAS/scripts/git_up/beta_software/summary_expr_beta.R $path1/Lung_h2_${h2}/${rep}/
Rscript /ysm-gpfs/pi/zhao-data/wl382/cWAS/scripts/git_up/beta_software/summary_beta_S.R $path1/Lung_h2_${h2}/${rep}/
#script summary_files.R /cWAS/data/simul/gwas/${rep}/
Rscript /ysm-gpfs/pi/zhao-data/wl382/cWAS/scripts/git_up/beta_software/summary_files.R $path1/Lung_h2_${h2}/${rep}/

# plink
plink --bfile /ysm-gpfs/pi/zhao-data/wl382/cWAS/data/GTEX/v6/geno/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs --extract $path1/Lung_h2_${h2}/${rep}/gtex_v8_all.bim --make-bed --out $path1/Lung_h2_${h2}/${rep}/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs_all

module load R/3.6.1-foss-2018b;
Rscript /ysm-gpfs/pi/zhao-data/wl382/cWAS/scripts/git_up/beta_software/assoc_test_v2.R --out=$path1/Lung_h2_${h2}/${rep}/assoc_LM22_v6 --signature-matrix=$path1/Lung_h2_${h2}/${rep}/LM22_new_all.txt --gwas-beta=$path1/Lung_h2_${h2}/${rep}/gwas_all.txt --ref-geno=$path1/Lung_h2_${h2}/${rep}/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs_all --ref-cellFraction=/ysm-gpfs/pi/zhao-data/wl382/cWAS/data/GTEX/v6/ref_cell_fracs_predExpr_Whole_Blood.txt --tol-inv=0.01 --betaS-file=$path1/Lung_h2_${h2}/${rep}/expr_beta_S_all.txt;
#rm $path1/Lung_h2_${h2}/${rep}/*bim;
#rm $path1/Lung_h2_${h2}/${rep}/*bed;
#rm $path1/Lung_h2_${h2}/${rep}/*fam;
#rm -rf $path1/Lung_h2_${h2}/${rep}/*chr*txt;
#rm -rf $path1/Lung_h2_${h2}/${rep}/*txt;
#cd $path1/Lung_h2_${h2}/${rep};
#rm -rf *log;rm -rf *nosex
