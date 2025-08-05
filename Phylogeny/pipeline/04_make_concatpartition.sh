#!/usr/bin/bash -l
#SBATCH -p short -c 2 --mem 16gb --out logs/make_mfa_pep.%A.log
CPU=${SLURM_CPUS_ON_NODE}
if [ -z $CPU ]; then
    CPU=1
fi
CPURUN=96

module load phykit
PREFIX=fungi5k_v1
USERTREE=fungi_tree/final_tree.nw
for type in fungi 
do
	FILTERDIR=${type}_msa_filter
	STEM=${type}
	if [ ! -d $FILTERDIR/selected_MSAs ]; then
		echo "no $FILTERDIR/selected_MSAs folder"
		continue
	fi
	pushd $FILTERDIR/selected_MSAs
	ls *.mfa > filenames
	mkdir -p ../../${FILTERDIR}-buildtree
	phykit create_concat -a filenames -p ../../$FILTERDIR-buildtree/${PREFIX}.${STEM}
	popd
	pushd $FILTERDIR-buildtree
	perl -i -p -e 's/AUTO/PROT/' ${PREFIX}.${STEM}.partition
	#sbatch -p stajichlab -c $CPURUN --mem 128gb -J modeltest$type --out modeltest-${type}.%A.log --wrap "module load modeltest-ng; modeltest-ng -i ${PREFIX}.${STEM}.fa -q ${PREFIX}.${STEM}.partition --processes $CPURUN -T raxml -d aa -t user -u $USERTREE -T raxml"
	sbatch -p stajichlab -c $CPURUN --mem 128gb -J modeltest$type --out modeltest-${type}.%A.log --wrap "hostname; module load modeltest-ng; modeltest-ng -i ${PREFIX}.${STEM}.fa -q ${PREFIX}.${STEM}.partition --processes $CPURUN -T raxml -d aa -t mp -T raxml"
	popd
done
