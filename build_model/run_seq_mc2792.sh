cp /ysm-gpfs/pi/zhao-data/wl382/cWAS/scripts/git_up/build_model/models_chr*_mc2792_remain.sh .
rm models_all_chrs_mc2792.sh
for i in {6..14};do
cat models_chr${i}_mc2792_remain.sh >> models_all_chrs_mc2792.sh
done
#cat models_chr*_mc2792.sh >> models_all_chrs_mc2792.sh
split --suffix-length=3 --numeric-suffixes --lines=50000  models_all_chrs_mc2792.sh models_all_chrs_mc2792_sh_
for i in {000..030}; do # change 030 to the corresponding largest number with the prefix models_all_chrs_zd72_sh_
  module load dSQ
	dSQ --jobfile models_all_chrs_mc2792_sh_${i} --mem-per-cpu=50g -t 1:00:00 -n 1 -p scavenge,general -J model
  sbatch dsq-models_all_chrs_mc2792_sh_${i}-2020-01-02.sh
done

