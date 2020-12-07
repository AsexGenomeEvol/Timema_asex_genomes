#!/bin/bash

#get variables to use in bash from two columns of a file per row

#filter RepeatMasker out.gff by identity percentage and length

#folders needed
#gff

#add modules

clear
while read line; do
  read -a twoSp <<< $line
  printf "\n****processing ${twoSp[0]}*****\n"
  SP1=${twoSp[0]}
  SP2=${twoSp[1]}


#filter
cd gff
#cat $SP1''_b3v07.fa.out.gff | awk '{if($6<=10)print($0);}' | awk '{if(($5-$4)>=80)print($0);}' > $SP1''_b3v07.fa.out.len80div10.gff
#cat $SP1''_b3v07.fa.out.gff | awk '{if($6<=5)print($0);}' | awk '{if(($5-$4)>=80)print($0);}' > $SP1''_b3v07.fa.out.len80div5.gff
#cat $SP1''_b3v07.fa.out.gff | awk '{if($6<=30)print($0);}' | awk '{if(($5-$4)>=200)print($0);}' > $SP1''_b3v07.fa.out.len200div30.gff
#cat $SP1''_b3v07.fa.out.gff | awk '{if($6<=10)print($0);}' | awk '{if(($5-$4)>=200)print($0);}' > $SP1''_b3v07.fa.out.len200div10.gff
#cat $SP1''_b3v07.fa.out.gff | awk '{if($6<=20)print($0);}' | awk '{if(($5-$4)>=80)print($0);}' > $SP1''_b3v07.fa.out.len80div20.gff
cat $SP1''_b3v08.fasta.out.gff | awk '{if(($5-$4)>=80)print($0);}' > $SP1''_b3v08.fa.out.len80.gff

cd ..

done < ./samplenames

printf "****DONE****\n"
