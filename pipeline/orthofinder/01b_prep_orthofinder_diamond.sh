#!/usr/bin/bash -l
#SBATCH -N 1 -c 16 --mem 24gb --out logs/orthofinder.%A.log

mkdir -p logs
module load orthofinder
opt="" # could change to "-C xeon" and will run on the xeon nodes; # could change this to empty and will run on any node
JOBS=orthofinder_steps.diamond.sh
LOG=orthofinder_steps.diamond.log
CHUNK=15000
export TEMPDIR=$SCRATCH
if [ ! -f $LOG ]; then
	orthofinder -op -t 16 -a 16 -f input -S diamond -o OrthoFinder_diamond > $LOG
fi
grep ^diamond $LOG | grep -v 'commands that must be run' | perl -p -e 's/-p 1/-p 8/g'> $JOBS
t=$(wc -l $JOBS | awk '{print $1}')
MAX=$(expr $t / $CHUNK)
echo "t is $t MAX is $MAX"
exit
for n in $(seq $MAX)
do
	START=$(perl -e "printf('%d',1 + $CHUNK * ($n - 1))")
	END=$(perl -e "printf('%d',$CHUNK* $n)")
#	echo "$START,$END for $n"
	run=$(sed -n ${START},${END}p $JOBS)
	sbatch -p short $opt --out logs/diamond.$n.log -J Dmd$n -N 1 -n 1 -c 8 --mem 4gb --wrap "module load orthofinder; $run"
done
