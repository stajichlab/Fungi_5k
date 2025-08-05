#!/usr/bin/bash -l
cd input
for a in $(ls -U); do
	m=$(head -n 1 $a | grep "^>" | perl -p -e 's/^>([^_]+)_.+/$1/'); 
	ln -s $(realpath $a) ../../short_names_pep/$m.fasta; 
done

