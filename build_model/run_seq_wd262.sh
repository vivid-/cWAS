cp /ysm-gpfs/pi/zhao-data/wl382/cWAS/scripts/git_up/build_model/models_chr*_wd262.sh .
cat models_chr*_wd262.sh >> models_all_chrs_wd262.sh
split --suffix-length=3 --numeric-suffixes --lines=50000  models_all_chrs_wd262.sh models_all_chrs_wd262_sh_
for i in {000..030}; do # change 030 to the corresponding largest number with the prefix models_all_chrs_wd262_sh_
  module load dSQ
	dSQ --jobfile models_all_chrs_wd262_sh_${i} --mem-per-cpu=50g -C avx2 -t 1:00:00 -n 1 -p scavenge,general -J model
  sbatch dsq-models_all_chrs_wd262_sh_${i}-2019-11-29.sh
done

