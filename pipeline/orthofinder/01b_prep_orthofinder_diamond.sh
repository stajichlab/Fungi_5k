#!/usr/bin/bash -l
#SBATCH -N 1 -c 16 --mem 24gb --out logs/orthofinder.%A.log

module load orthofinder
mkdir -p logs
opt="" # could change to "-C xeon" and will run on the xeon nodes; # could change this to empty and will run on any node
JOBS=orthofinder_steps.diamond.sh
LOG=orthofinder_steps.diamond.log
CHUNK=100
export TEMPDIR=$SCRATCH
INDIR=input
if [ ! -f $LOG ]; then
	orthofinder -op -t 16 -a 16 -f $INDIR -S diamond -o OrthoFinder_diamond > $LOG
fi
grep ^diamond $LOG | grep -v 'commands that must be run' | perl -p -e 's/-p 1/-p 8/g'> $JOBS
t=$(wc -l $JOBS | awk '{print $1}')
MAX=$(expr $t / $CHUNK + 1) # + 1 to round up
echo "t is $t MAX is $MAX"

# COULD FIGURE OUT HOW TO DO THIS AS ARRAY JOB FOR EASIER START/DELETE JOBS
for n in $(seq $MAX)
do
	START=$(perl -e "printf('%d',1 + $CHUNK * ($n - 1))")
	END=$(perl -e "printf('%d',$CHUNK * $n)")
	if [[ $START -gt $t ]]; then
		START=$t
	fi
	if [[ $END -gt $t ]]; then
		END=$t
	fi
	#END=$(expr $START + 1) # debugging set this to + 1
	echo "$START,$END for $n"
# this can be uncommented for some reporting and checking on what needs to be rerun
#	for line in $(seq $START 1 $END)
#	do 
#		runline=$(sed -n ${line}p $JOBS)
#		fileis=$(sed -n ${line}p $JOBS | awk '{print $9}')
#		if [ ! -f $fileis.gz ]; then
#			echo "need to run $fileis -> Running $runline"
#		fi
#		echo ${line}
#	done

	# this supports re-running jobs without redoing all that are already finished, 
	# old loop in the wrap code was 
	# for line in \$(seq $START 1 $END); do time \$(sed -n \${line}p $JOBS); echo ${line}; date; done"
	sbatch $opt --out logs/diamond.$n.log -J Dmd$n -N 1 -n 1 -c 12 --mem 4gb --wrap "module load orthofinder; \
		for line in \$(seq $START 1 $END); do fileis=\$(sed -n \${line}p $JOBS | awk '{print \$9}'); if [ ! -f \$fileis.gz ]; then time \$(sed -n \${line}p $JOBS); echo ${line}; date; fi; done"
done
