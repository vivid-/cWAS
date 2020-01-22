#051
for i in {016..030}; do
  module load dSQ
	dSQ --jobfile models_sh_${i} --mem-per-cpu=50g -t 1:00:00 -n 1 -p scavenge,general -J model
  sbatch dsq-models_sh_${i}-2019-11-27.sh
done
