#!/bin/bash

#get variables to use in bash from two columns of a file per row

#map reads to timema genomes

#folder needed:
# gff
# map
# count


#add modules
module add UHTS/Analysis/HTSeq/0.9.1;


clear
while read line; do
  read -a twoSp <<< $line
  printf "\n****processing ${twoSp[0]}*****\n"
  SP1=${twoSp[0]}


#HTSeq-count
printf "**HTseq-count $SP1**\n\n"
#htseq-count -f bam -r name -s no -t similarity -i Target --nonunique none map/$SP1''_aln.sorted.bam gff/$SP1''_b3v06.fa.out.len80div5.gff &> ./count_len80div5/$SP1''_htseq_unq.cnt
htseq-count -f bam -r name -s no -t similarity -i Target --nonunique none map/$SP1''_aln.sorted.bam gff/$SP1''_b3v08.fa.out.len80.gff &> ./count/$SP1''_htseq_unq.cnt


done < ./samplenames

module rm UHTS/Analysis/HTSeq/0.9.1;


printf "****DONE****\n"



#GETTING READS HIT PER ELEMENT
#use RM gff converter on masked.out files
#rmOutToGFF3.pl
# /software/SequenceAnalysis/Repeat/RepeatMasker/4.0.7/util/rmOutToGFF3.pl

#for a better summary table:
# /software/SequenceAnalysis/Repeat/RepeatMasker/4.0.7/util/buildSummary.pl
