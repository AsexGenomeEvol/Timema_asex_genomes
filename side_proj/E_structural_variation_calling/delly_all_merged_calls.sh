#!/bin/bash

# this is both merging AND filtering
SP=$1

delly merge -o data/$SP/variant_calls/"$SP"_all_calls_merged.bcf \
    data/$SP/variant_calls/*/*_delly.bcf \
    data/$SP/variant_calls/*/*_smoove/*-smoove.genotyped.vcf \
    data/$SP/variant_calls/*/*_manta/results/variants/diploidSV.vcf
# takes less then minute
