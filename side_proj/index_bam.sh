#!/bin/bash
#
# index bam specified by 1st argument

bsub <<< """
#BSUB -L /bin/bash
#BSUB -J index_bam
#BSUB -q priority
#BSUB -o index_$1.log
#BSUB -n 1
#BSUB -M 1000000

module add UHTS/Analysis/samtools/1.3

samtools index $1
"""
