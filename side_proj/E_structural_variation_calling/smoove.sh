#!/bin/bash

# E_structural_variation_calling/smoove.sh
#   data/1_Tdi/reference/1_Tdi_b3v08.fasta.gz
#   data/1_Tdi/reference/1_Tdi_b3v08.fasta.gz.fai
#   data/1_Tdi/reference/1_Tdi_b3v08.fasta.gz.gzi
#   data/1_Tdi/mapping/Tdi_01_to_b3v08_mapped_within_scfs.bam
#   data/1_Tdi/mapping/Tdi_01_to_b3v08_mapped_within_scfs.bam.bai
#   data/1_Tdi/variant_calls/Tdi_01/Tdi_01_smoove
SAMPLE=$1
FASTA=$2
BAM=$5
OUTDIR=$7

source /scratch/beegfs/monthly/kjaron/src/svtyper/bin/activate
module add UHTS/Analysis/samtools/1.8

smoove call -x --genotype --name $SAMPLE -f $FASTA --processes 8 $BAM --outdir $OUTDIR
