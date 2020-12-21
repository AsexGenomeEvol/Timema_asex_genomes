#!/bin/bash

for sp in $TIMEMAS; do
    for ind in 00 01 02 03 04 05; do
        DIR=data/$sp/trimmed_reads/$(echo $sp | cut -f 2 -d _)_"$ind"/
        bsub -o data/$sp/coverage_"$ind".log 'find '$DIR' -name "*q.gz" -exec scripts/generic_genomics/fastq.gz2number_of_nt.sh {} \; > '"$DIR""$sp"_coverage.stats
    done
done