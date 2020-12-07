#!/bin/bash
module add UHTS/Analysis/usearch/10.0.240;


clear
while read line; do
  read -a twoSp <<< $line
  printf "\n****processing ${twoSp[0]}*****\n"
  SP1=${twoSp[0]}


 printf "**running $SP1**\n\n"


cat $SP1''_all_repeats.fasta.rfmt $SP1''.consensi.fa.classified > $SP1''_merged.TElib.fa

usearch -cluster_fast $SP1''_merged.TElib.fa -id 0.95 -centroids $SP1''_merged.TElib.fa.centroid95 -uc $SP1''_merged.TElib.fa.clusters95.uc -threads 24
#usearch -cluster_fast $SP1''_merged.TElib.fa -id 0.80 -centroids $SP1''_merged.TElib.fa.centroid80 -uc $SP1''_merged.TElib.fa.clusters80.uc -threads 24

done < samplenames

module rm UHTS/Analysis/usearch/10.0.240;

printf "****DONE****\n"

