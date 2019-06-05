#!/bin/bash
SP=$1

GENOME=data/$SP/reference/"$SP"_b3v08.fasta.gz
MERGED_BCF=data/$SP/variant_calls/all_calls_merged.bcf

for BAM in data/$SP/mapping/*_mapped_within_scfs.bam; do
    SAMPLE=$(echo $BAM | cut -f 4 -d / | cut -f 1,2 -d _)
    OUTPUT=data/"$SP"/variant_calls/"$SAMPLE"/delly_genotyped_error_candidates.bcf
    delly call -g $GENOME -v $MERGED_BCF -o $OUTPUT $BAM &
done

wait
echo $SP "Done!"


