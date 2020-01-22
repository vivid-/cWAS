for i in {1..5};do
  cp /ysm-gpfs/pi/zhao-data/wl382/cWAS/scripts/git_up/build_model/models_chr${i}_bz234.sh .
done
rm models_all_chrs_bz234.sh
cat models_chr*_bz234.sh >> models_all_chrs_bz234.sh
split --suffix-length=3 --numeric-suffixes --lines=50000  models_all_chrs_bz234.sh models_all_chrs_bz234_sh_
for i in {000..008}; do # change 030 to the corresponding largest number with the prefix models_all_chrs_zd72_sh_
  module load dSQ
	dSQ --jobfile models_all_chrs_bz234_sh_${i} --mem-per-cpu=50g -t 1:00:00 -n 1 -p scavenge,general -J model
  sbatch dsq-models_all_chrs_bz234_sh_${i}-2019-12-12.sh
done

