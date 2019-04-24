#!/bin/bash
#
# SP
# IS

# crashed because of memory on local (now it computes directly to beegfs)
# it was about to crush because of memory, but it was maybe because it knew that it has exactly 20G

#BSUB -L /bin/bash
#BSUB -J $1_dedup
#BSUB -q normal
#BSUB -o dedupl.out
#BSUB -e dedupl.err
#BSUB -n 1
#BSUB -M 40000000

module add UHTS/Analysis/picard-tools/2.2.1

# make local directory for computations
LOCAL_DIR=/scratch/local/monthly/kjaron/$1/$2
export TMP_DIR=\$LOCAL_DIR/temp
mkdir -p \$LOCAL_DIR/temp

# run sniffles
picard-tools MarkDuplicates \
  INPUT=$TROOT/data/$1/variant_calling/GATK/$2/map_pe_to_$1.bam \
  OUTPUT=dedup_pe_to_$1.bam \
  METRICS_FILE=metrics.txt
