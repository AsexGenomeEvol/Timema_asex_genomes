#!/bin/bash

#map reads to timema genomes v8

#folder needed:
# reads
# map
# genomes


########
#BEFORE#
########


clear
while read line; do
  read -a twoSp <<< $line
  printf "\n****processing ${twoSp[0]}*****\n"
  SP1=${twoSp[0]}
  SP2=${twoSp[1]}



#HTSeq-count
printf "**HTseq-count $SP1**\n\n"
htseq-count -f bam -r name -s no -t similarity -i Target --nonunique none map/$SP1$SP2''_aln.sorted.bam gff/$SP1''_b3v08.fa.out.len80.gff &> ./count/$SP1$SP2''_htseq_unq.cnt


done < ./filelist.htseq2


printf "****DONE****\n"
