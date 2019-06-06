#!/bin/bash
SP=$1

GENOME=data/$SP/reference/"$SP"_b3v08.fasta.gz
MERGED_BCF=data/$SP/variant_calls/"$SP"_survivor_all_calls_union.vcf.gz

bcftools index $MERGED_BCF

for BAM in data/$SP/mapping/*_mapped_within_scfs.bam; do
    SAMPLE=$(echo $BAM | cut -f 4 -d / | cut -f 1,2 -d _)
    OUTPUT=data/"$SP"/variant_calls/"$SAMPLE"/"$SAMPLE"_delly_genotyping.bcf
    delly call -g $GENOME -v $MERGED_BCF -o $OUTPUT $BAM &
done

wait

bcftools merge -m id -O b -o data/"$SP"/variant_calls/"$SP"_delly_genotyping_merged.bcf data/"$SP"/variant_calls/*/delly_genotyping.bcf
bcftools view data/"$SP"/variant_calls/"$SP"_delly_genotyping_merged.bcf > data/"$SP"/variant_calls/"$SP"_delly_genotyping_merged.vcf

echo $SP "Done!"


