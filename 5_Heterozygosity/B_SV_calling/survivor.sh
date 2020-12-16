#!/bin/bash

for SP in 4_Tte 4_Tbi 3_Tms 3_Tce 1_Tdi 1_Tps 2_Tsi 2_Tcm 5_Tge 5_Tpa; do
    echo $SP
    BASE=data/manta_SV_calls/data/$SP

    ls $BASE/*_manta/results/variants/diploidSV_reduced.vcf > $BASE/variant_calls.tsv
    ls $BASE/*_manta/results/variants/diploidSV_filt_relaxed.vcf > $BASE/variant_calls_filt_relaxed.tsv
    ls $BASE/*_manta/results/variants/diploidSV_filt_stringent.vcf > $BASE/variant_calls_filt_stringent.tsv
    ls $BASE/*_manta/results/variants/diploidSV_filt_very_stringent.vcf > $BASE/variant_calls_filt_very_stringent.tsv

    # SURVIVOR stats "$MANTA" 1 500000 5 "$BASE"_manta.summary
    # SURVIVOR stats "$LUMPY" 1 500000 5 "$BASE"_lumpy.summary
    # SURVIVOR stats "$DELLY" 1 500000 5 "$BASE"_delly.summary

    SURVIVOR merge $BASE/variant_calls.tsv 100 1 1 1 0 30 "$BASE"/SVs_no_filt_union.vcf
    SURVIVOR merge $BASE/variant_calls_filt_relaxed.tsv 100 1 1 1 0 30 "$BASE"/SVs_filt_relaxed_union.vcf
    SURVIVOR merge $BASE/variant_calls_filt_stringent.tsv 100 1 1 1 0 30 "$BASE"/SVs_filt_stringent_union.vcf
    SURVIVOR merge $BASE/variant_calls_filt_very_stringent.tsv 100 1 1 1 0 30 "$BASE"/SVs_filt_very_stringent_union.vcf


    SURVIVOR stats "$BASE"/SVs_no_filt_union.vcf -1 -1 -1 "$BASE"/SVs_no_filt_union.stats
    SURVIVOR stats "$BASE"/SVs_filt_relaxed_union.vcf -1 -1 -1 "$BASE"/SVs_filt_relaxed_union.stats
    SURVIVOR stats "$BASE"/SVs_filt_stringent_union.vcf -1 -1 -1 "$BASE"/SVs_filt_stringent_union.stats
    SURVIVOR stats "$BASE"/SVs_filt_very_stringent_union.vcf -1 -1 -1 "$BASE"/SVs_filt_very_stringent_union.stats
done