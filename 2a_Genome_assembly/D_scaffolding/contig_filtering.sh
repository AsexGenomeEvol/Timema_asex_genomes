#!/bin/bash
#
# the first argument is the assembly

if [[ ! -s $1 ]]
then
	>&2 echo GENOME $1 IS MISSING;
	exit 1
fi

bsub <<< """
#BSUB -L /bin/bash
#BSUB -J $(basename $2)
#BSUB -q priority
#BSUB -n 1
#BSUB -M 1000000

python3 $TROOT/scripts/generic_genomics/fasta2fasta_length_filtering.py $1 250 > $2
"""
