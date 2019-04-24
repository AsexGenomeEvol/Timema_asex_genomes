#!/bin/bash
#
# build dict from fasta specified by 1st argument

#BSUB -L /bin/bash
#BSUB -J dict_fa
#BSUB -q normal
#BSUB -o dict_fa.out
#BSUB -e dict_fa.err
#BSUB -n 1
#BSUB -M 4000000

module add UHTS/Analysis/picard-tools/2.2.1

FILENAME=$1
picard-tools CreateSequenceDictionary R=$1 O="\${FILENAME%.*}".dict
