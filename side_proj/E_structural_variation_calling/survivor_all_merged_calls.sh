#!/bin/bash

module add UHTS/Analysis/SURVIVOR/1.0.2

SP=$1

# data/$SP/variant_calls/*/*_delly.bcf \
ls data/$SP/variant_calls/*/*_smoove/*-smoove.genotyped.vcf \
   data/$SP/variant_calls/*/*_manta/results/variants/diploidSV.vcf \
   > data/$SP/variant_calls/all_sv_files

SURVIVOR merge data/$SP/variant_calls/all_sv_files 1000 1 1 1 0 30 data/$SP/variant_calls/"$SP"_survivor_all_calls_union.vcf


# I used [SURVIVOR][1] to merge structural variant (SV) calls from several different SV callers and individuals. Now I have a catalogue of candidate that I would like to genotype in the same set of individuals. It's basically, what authors of Delly propose, but expanded to multiple SV callers.
#
# ```
# SURVIVOR merge all_sv_files 1000 1 1 1 0 30 sv_calls_union
# ```
#
#   [1]: https://github.com/fritzsedlazeck/SURVIVOR