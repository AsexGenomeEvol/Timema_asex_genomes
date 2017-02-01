#!/bin/bash
#
# species
# insert size
# ref name

bsub <<< """
#BSUB -L /bin/bash
#BSUB -J callBayes_$1
#BSUB -q priority
#BSUB -o atlas.out
#BSUB -e index_fa.err
#BSUB -n 1
#BSUB -M 1000000

module add UHTS/Analysis/samtools/1.3

# make local directory for computations
LOCAL_DIR=/scratch/local/monthly/kjaron/$1/$2
export TMP_DIR=\$LOCAL_DIR/temp
mkdir -p \$LOCAL_DIR/temp

atlas task=callBayes \
  bam=\$TROOT/data/$1/variant_calling/GATK/$2/map_pe_to_$1.bam \
  1> callBayes.out 2> callBayes.err¬
"""
