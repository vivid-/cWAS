#step1: get the weights of those signature genes
Rscript signaturea_weights_v8.R $sig_file $tissue
#step2: impute gene expression based for those genes
Rscript generate_bash_pred_expr_sh.R $sig_file $tissue
module load dSQ
dSQ --jobfile pred_expr_v8_${tissue}.sh --mem-per-cpu=100g -t 4:00:00 -n 1 -p scavenge,general,pi_zhao -J wht

#step3: summarize the predicted gene expression
Rscript summary_expr_pred_v8.R $sig_file $tissue
#step4: estimate the cell fractions based on predicted gene expression
Rscript linear_fracs.R $sig_file $tissue
