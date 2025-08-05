#!/usr/bin/bash -l
#SBATCH -p epyc -c 40 --mem 128gb --out logs/fungi5k_raxml.log

module load raxml-ng

pushd fungi_msa_filter-buildtree

raxml-ng --all --msa fungi5k_v1.fungi.fa.raxml.rba --tree pars{10} -d aa --bs-trees 200 --threads auto{40} --workers auto{4}
