#!/bin/bash
#make repeat landscapes from RepeatMasker output


clear
while read line; do
  read -a twoSp <<< $line
  printf "\n****processing ${twoSp[0]}*****\n"
  SP1=${twoSp[0]}
  SP2=${twoSp[1]}


  printf "**running species $SP1 with genome size $SP2**\n\n"

#masking

cd $SP1''_RM_div30;
perl /NVME/Software/Repeats/RepeatMasker/util/calcDivergenceFromAlign.pl -s $SP1''_b3v08.fa.divsum $SP1''_b3v08.fasta.align
perl /home/jbast/Scratch/jbast/timema/masking/createRepeatLandscape_mod.pl -div $SP1''_b3v08.fa.divsum -g $SP2 > $SP1''_b3v08.fa_mod.html
cd ..

#Patrick script
cd $SP1''_RM_div30;
python /home/jbast/Scratch/jbast/timema/masking/te.py -s repeatmasker2r -i1 $SP1''_b3v08.fa_mod.html -o $SP1''_b3v08_R_data.txt
cd ..


done < samplenames

printf "****DONE****\n"


