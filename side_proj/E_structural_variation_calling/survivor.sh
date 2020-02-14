#!/bin/bash

#module add UHTS/Analysis/SURVIVOR/1.0.2
#use THIS instead

/scratch/temporary/ptranvan/Software/survivor/1.0.6/Debug/SURVIVOR -h

# gather of all vcf files
# sort them??

SP=3_Tms
# by SAMPLE
# SAMPLE=Tms_01

for SAMPLE in $(ls variant_calls); do
    echo $SAMPLE

    DELLY="variant_calls/"$SAMPLE"/"$SAMPLE"_delly.vcf"
    MANTA="variant_calls/"$SAMPLE"/"$SAMPLE"_manta/results/variants/diploidSV.vcf"
    LUMPY="variant_calls/"$SAMPLE"/"$SAMPLE"_smoove/"$SAMPLE"-smoove.genotyped.vcf"

    # DELLY="variant_calls/"$SAMPLE"/"$SAMPLE"_delly.bcf"
    # MANTA="variant_calls/"$SAMPLE"/"$SAMPLE"_manta/results/variants/diploidSV.vcf.gz"
    # LUMPY="variant_calls/"$SAMPLE"/"$SAMPLE"_smoove/"$SAMPLE"-smoove.genotyped.vcf.gz"
    # gunzip $MANTA
    # gunzip $LUMPY
    # MANTA=${MANTA%.*}
    # LUMPY=${LUMPY%.*}
    # bcftools view $DELLY > ${DELLY%.*}.vcf
    # DELLY=${DELLY%.*}.vcf
    BASE=variant_calls/"$SAMPLE"/"$SAMPLE"

    SURVIVOR stats "$MANTA" 1 500000 5 "$BASE"_manta.summary
    SURVIVOR stats "$LUMPY" 1 500000 5 "$BASE"_lumpy.summary
    SURVIVOR stats "$DELLY" 1 500000 5 "$BASE"_delly.summary

    ls $DELLY $MANTA $LUMPY > sample_files

    SURVIVOR merge sample_files 1000 1 1 1 0 30 "$BASE"_union.vcf
    SURVIVOR merge sample_files 1000 2 1 1 0 30 "$BASE"_merged.vcf
    SURVIVOR merge sample_files 1000 3 1 1 0 30 "$BASE"_strict.vcf

    SURVIVOR stats "$BASE"_union.vcf -1 -1 -1 "$BASE"_union.summary
    SURVIVOR stats "$BASE"_merged.vcf -1 -1 -1 "$BASE"_merged.summary
    SURVIVOR stats "$BASE"_strict.vcf -1 -1 -1 "$BASE"_strict.summary
done

## merge by soft
ls variant_calls/*/*_delly.vcf > delly_files
ls variant_calls/*/*_manta/results/variants/diploidSV.vcf > manta_files
ls variant_calls/*/*_smoove/*-smoove.genotyped.vcf > lumpy_files

for soft in delly manta lumpy; do
    BASE=variant_calls/"$soft"/merged
    mkdir -p variant_calls/"$soft"

    SURVIVOR merge "$soft"_files 1000 1 1 1 0 30 "$BASE"_union.vcf
    SURVIVOR merge "$soft"_files 1000 2 1 1 0 30 "$BASE"_merged.vcf
    SURVIVOR merge "$soft"_files 1000 6 1 1 0 30 "$BASE"_all.vcf

    SURVIVOR stats "$BASE"_union.vcf -1 -1 -1 "$BASE"_union.summary
    SURVIVOR stats "$BASE"_merged.vcf -1 -1 -1 "$BASE"_merged.summary
    SURVIVOR stats "$BASE"_all.vcf -1 -1 -1 "$BASE"_all.summary
done

ls variant_calls/*/*_all.vcf > sample_files

SURVIVOR merge sample_files 1000 1 1 1 0 30 asm_errors_union.vcf
SURVIVOR merge sample_files 1000 2 1 1 0 30 asm_errors_merged.vcf
SURVIVOR merge sample_files 1000 3 1 1 0 30 asm_errors_strict.vcf

SURVIVOR stats asm_errors_union.vcf -1 -1 -1 asm_errors_union.summary
SURVIVOR stats asm_errors_merged.vcf -1 -1 -1 asm_errors_merged.summary
SURVIVOR stats asm_errors_strict.vcf -1 -1 -1 asm_errors_strict.summary