#!/bin/bash

for ASSEMBLY in "$@"
do
	echo "$ASSEMBLY"
	python $TROOT/scripts/generic_genomics/fasta2genomic_stats.py $ASSEMBLY 1300000000
done
