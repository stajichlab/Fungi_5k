#!/usr/bin/bash -l
mkdir -p cds pep GFF

for a in $(grep Sordariomycetes ncbi_accessions_taxonomy.csv  | cut -d, -f1); do 
	in=$(ls /bigdata/stajichlab/shared/projects/1KFG/2021/source/NCBI/pep/${a}*)
	if [ ! -z "$in" ]; then ln -s $in pep/ ; fi
	in=$(ls /bigdata/stajichlab/shared/projects/1KFG/2021/source/NCBI/cds/${a}*)
	if [ ! -z "$in" ]; then ln -s $in cds/ ; fi
	in=$(ls /bigdata/stajichlab/shared/projects/1KFG/2021/source/NCBI/GFF/${a}*)
	if [ ! -z "$in" ]; then 
		ln -s $in GFF/ ; 
	 	ln -s /bigdata/stajichlab/shared/projects/1KFG/2021/source/NCBI/DNA/${a}*  DNA/
	fi
done

