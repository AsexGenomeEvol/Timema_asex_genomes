#!/bin/bash
#
# index fasta specified by 1st argument

bsub <<< """
#BSUB -L /bin/bash
#BSUB -J index_fa
#BSUB -q priority
#BSUB -o index_fa.out
#BSUB -e index_fa.err
#BSUB -n 1
#BSUB -M 1000000

module add UHTS/Analysis/samtools/1.3

samtools faidx $1
"""
