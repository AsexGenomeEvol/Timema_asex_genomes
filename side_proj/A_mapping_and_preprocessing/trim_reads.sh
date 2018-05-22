#!/bin/bash

# This is a script that will trimm one by one all the reads of one re-seq individual
# SP
# SAM
# the other arguments should be
# - directory with raw reads
# - the output directory with trimmed reads

module add UHTS/Analysis/trimmomatic/0.36;

SP=$1
SAM=$2
ADAPTERS=$3
# data/4_Tbi/raw_reads/Tbi_01/
INPUT_DIR=$4
# data/4_Tbi/trimmed/Tbi_01/
OUTPUT_DIR=$5

R1_READS=$(ls $INPUT_DIR | grep "_R1_")

mkdir $OUTPUT_DIR

for r1 in $R1_READS; do
    r2=${r1/R1/R2} # get R2 reads by substituting R1 in the r1 name
    NAMEBASE=$(basename $r1 .fastq.gz | awk -F "_R1" '{print $1 $2}') # delete suffix and _R? from the read names
    t_r1="$NAMEBASE"_R1.cleaned.fastq.gz # specify output name based on the input
    t_r2="$NAMEBASE"_R2.cleaned.fastq.gz # naminch scheme taken from Emeric
    np_r1="$NAMEBASE"_R1.se.cleaned.fastq.gz
    np_r2="$NAMEBASE"_R2.se.cleaned.fastq.gz

    trimmomatic PE -threads 16 "$INPUT_DIR"/"$r1" "$INPUT_DIR"/"$r2" \
        "$OUTPUT_DIR"/$t_r1 "$OUTPUT_DIR"/$np_r1 \
        "$OUTPUT_DIR"/$t_r2 "$OUTPUT_DIR"/$np_r1 \
        ILLUMINACLIP:"$ADAPTERS":3:25:6 LEADING:9 TRAILING:9 SLIDINGWINDOW:4:15 MINLEN:96
done
