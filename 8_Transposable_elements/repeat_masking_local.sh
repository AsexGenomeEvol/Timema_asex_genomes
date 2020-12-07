#!/bin/bash
#RepeatMasker run on Motoko server

clear
while read line; do
  read -a twoSp <<< $line
  printf "\n****processing ${twoSp[0]}*****\n"
  SP1=${twoSp[0]}



  printf "**running $SP1**\n\n"

#masking
cd /home/jbast/Scratch/jbast/timema/masking;
mkdir $SP1''_RM_div30;


/home/jbast/Software/Repeats/RepeatMasker/RepeatMasker -lib /home/jbast/Scratch/jbast/timema/repeat_libs/Timema_TElibs.merged.fa \
-gccalc -gff -u -a -xsmall -no_is -div 30 -engine rmblast -pa 60 \
/home/jbast/Scratch/jbast/timema/masking/genomes_v08/$SP1''_b3v08.fasta \
-dir /home/jbast/Scratch/jbast/timema/masking/$SP1''_RM_div30


done < samplenames

printf "****DONE****\n"
