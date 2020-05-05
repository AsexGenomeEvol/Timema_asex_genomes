#!/bin/bash

conda activate default_genomics

SCRIPT=scripts/convertManta/convertManta2Paragraph_compatible_vcf.py

for SP in 4_Tte 4_Tbi 3_Tms 3_Tce 1_Tdi 1_Tps 2_Tsi 2_Tcm 5_Tge 5_Tpa; do
    SP_NAME=$(echo $SP | cut -c3-5)
    REF=data/final_references/"$SP"_b3v08.fasta
    for ind in 00 01 02 03 04 05; do
        echo $SP $ind
        VARIANTS=data/manta_SV_calls/"$SP_NAME"_"$ind"_manta/results/variants/diploidSV.vcf.gz
        OUT=data/manta_SV_calls/"$SP_NAME"_"$ind"_manta/results/variants/diploidSV_paragraph_compatible.vcf
        python3 $SCRIPT $REF $VARIANTS > $OUT
    done
done
