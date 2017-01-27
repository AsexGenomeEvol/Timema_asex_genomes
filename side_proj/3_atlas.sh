#!/bin/bash
#
# index fasta specified by 1st argument

bsub <<< """
#BSUB -L /bin/bash
#BSUB -J index_fa
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

GenomeAnalysisTK \
  -T HaplotypeCaller \
  -R \$TROOT/data/$1/reference/$3 \
  -I \$TROOT/data/$1/variant_calling/GATK/$2/map_pe_to_$1.bam \
  --genotyping_mode DISCOVERY \
  -stand_call_conf 30 \
  -o raw_variants.vcf
"""
