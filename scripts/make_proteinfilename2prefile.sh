#!/usr/bin/bash -l
pushd input
grep -m 1 "^>" *.fa | perl -p -e 's/\.proteins\.fa//; s/:/\t/;s/>([^_]+)_.+/$1/' > ../proteinfile2prefix.tsv
