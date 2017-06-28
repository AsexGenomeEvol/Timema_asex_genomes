#!/bin/bash
#
# index fasta specified by 1st argument

bsub <<< """
#BSUB -L /bin/bash
#BSUB -J bwa_index_fa
#BSUB -q bgee
#BSUB -o bwa_ix_fa.log
#BSUB -e bwa_ix_fa.err
#BSUB -n 1
#BSUB -M 20000000

module add UHTS/Aligner/bwa/0.7.13

bwa index $1
"""
