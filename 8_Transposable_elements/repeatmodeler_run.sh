#!/bin/bash
#module add SequenceAnalysis/Repeat/RepeatMasker/4.0.5;
module add SequenceAnalysis/Repeat/RepeatModeler/1.0.8;



clear
while read line; do
  read -a twoSp <<< $line
  printf "\n****processing ${twoSp[0]}*****\n"
  SP1=${twoSp[0]}


 printf "**running $SP1**\n\n"


#compute on deeserv04 locally speeds up a lot

#cd /scratch/local/jbast/;
mkdir $SP1''_RepMod;

#cp /scratch/beegfs/monthly/kjaron/timema_assembly/data/$SP1/reference/$SP1''_b3v06.fa .



cd $SP1''_RepMod;

BuildDatabase -name $SP1 -engine ncbi ../../genomes/$SP1''_b3v06.fa

RepeatModeler -engine ncbi -pa 24 -database $SP1

cd ..

done < samplenames_test


module rm SequenceAnalysis/Repeat/RepeatModeler/1.0.8;
#module rm SequenceAnalysis/Repeat/RepeatMasker/4.0.5;

printf "****DONE****\n"

