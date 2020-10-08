#!/bin/bash
#
# SP
# IS
# ref name

#BSUB -L /bin/bash
#BSUB -J $1_GATK
#BSUB -q dee
#BSUB -o GATK.out
#BSUB -e GATK.err
#BSUB -n 1
#BSUB -M 20000000

module add UHTS/Analysis/GenomeAnalysisTK/3.7

# make local directory for computations
LOCAL_DIR=/scratch/local/monthly/kjaron/$1/$2
export TMP_DIR=\$LOCAL_DIR/temp
mkdir -p \$LOCAL_DIR/temp

GenomeAnalysisTK \
  -T HaplotypeCaller \
  -R \$TROOT/data/$1/reference/$3 \
  -I \$TROOT/data/$1/variant_calling/GATK/$2/map_pe_to_$1.bam \
  --genotyping_mode DISCOVERY \
  -stand_call_conf 30 \
  -o raw_variants.vcf
