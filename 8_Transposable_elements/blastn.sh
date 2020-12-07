#!/bin/bash

#blast TE libraries to Timema genomes

#folder needed:
# repeat_libraries_mod
# genomes


#add modules
module add Blast/ncbi-blast/2.7.1+;



clear
while read line; do
  read -a twoSp <<< $line
  printf "\n****processing ${twoSp[0]}*****\n"
  SP1=${twoSp[0]}

#blast prep
#for f in *.fa; do makeblastdb -dbtype nucl -in $f; done


#blastn
blastn -db ../genomes/$SP1''_b3v06.fa -query ../repeat_libraries_mod/$SP1''_merged.TElib.fa.centroid95min500.annot.srt.nbr.mod \
-task blastn -num_threads 24 -outfmt 6 -max_target_seqs 1 -out $SP1''.rep.blastn

#sort and merge
cat $SP1''.rep.blastn | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > $SP1''.rep.blastn2



done < ./samplenames

module rm Blast/ncbi-blast/2.7.1+;

printf "****DONE****\n"


