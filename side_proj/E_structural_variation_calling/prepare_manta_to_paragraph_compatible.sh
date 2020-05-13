#!/bin/bash

SCRIPT=scripts/convertManta/convertManta2Paragraph_compatible_vcf.py

rm E_structural_variation_calling/paragraph_commands.txt

for SP in 4_Tte 4_Tbi 3_Tms 3_Tce 1_Tdi 1_Tps 2_Tsi 2_Tcm 5_Tge 5_Tpa; do
# for SP in 4_Tte 3_Tms 1_Tdi 2_Tsi 5_Tge; do # this is just for asexuals (I want to staart with genotyping those)
# for SP in 4_Tbi 3_Tce 1_Tps 2_Tcm 5_Tpa; do # this is just for asexuals (I want to staart with genotyping those)
    SP_NAME=$(echo $SP | cut -c3-5)
    REF=data/final_references/"$SP"_b3v08.fasta
    for ind in 00 01 02 03 04 05; do
        echo $SP $ind
        VARIANTS=data/manta_SV_calls/"$SP_NAME"_"$ind"_manta/results/variants/diploidSV.vcf.gz
        OUT=data/manta_SV_calls/"$SP_NAME"_"$ind"_manta/results/variants/diploidSV_paragraph_compatible.vcf
        # prepare the vcf file
        python3 $SCRIPT $REF $VARIANTS > $OUT
        # prepare the paragraph genotyping command
        echo "multigrmpy.py -i $OUT -m data/genotyping/"$SP"_samples.txt -r $REF -o data/genotyping/"$SP"_ind"$ind"_genotyping --scratch-dir /scratch/kjaron/"$SP"_ind"$ind"_genotyping_temp --threads 8 && rm -r /scratch/kjaron/"$SP"_ind"$ind"_genotyping_temp" >> E_structural_variation_calling/paragraph_commands.txt
    done
done

