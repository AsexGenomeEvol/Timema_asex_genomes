#!/bin/bash

#get variables to use in bash from two columns of a file per row

#map reads to timema genomes

#folder needed:
# reads
# count
# map
# genomes


#get number of mapped reads from bam file to make input
#cd map
#for f in *.sorted.bam; do samtools view -c -F 260 $f; done > filestats.map
#cd ..
#paste samplenames map/filestats.map > samplenames_map_stats


clear
while read line; do
  read -a twoSp <<< $line
  printf "\n****processing ${twoSp[0]}*****\n"
  SP1=${twoSp[0]}
  SP2=${twoSp[1]}

#change folder CHANGE below as well
cd count;

#count
cat $SP1''_htseq_unq.cnt | grep 'Motif' | grep -v '(' | grep -v 'A-rich' | grep -v 'G-rich' | grep -v 'polypurine' | sed 's/"//g' | sed 's/Motif://g' | sed 's/ /\t/g' | cut -f1,4 | \
gawk '{sum[$1] += $2; N[$1]++ } END { for (key in sum) {summed = sum[key];printf "%s %d\n", key, summed; } }' | sort | sed 's/ /\t/g' | \
awk -v var="$SP2"  '{print $0 "\t" var}' | awk -v sp="$SP1" '{printf("%s\t%.6f\t%s\n",$0,$2/$3,sp)}' > $SP1''.unq.count

#cat $SP1''_htseq_nonunq.cnt | grep 'Motif' | grep -v '(' | grep -v 'A-rich' | grep -v 'G-rich'| grep -v 'polypurine' | sed 's/"//g' | sed 's/Motif://g' | sed 's/ /\t/g' | cut -f1,4 | \
#gawk '{sum[$1] += $2; N[$1]++ } END { for (key in sum) {summed = sum[key];printf "%s %d\n", key, summed; } }' | sort | sed 's/ /\t/g' | \
#awk -v var="$SP2"  '{print $0 "\t" var}' | awk -v sp="$SP1" '{printf("%s\t%.6f\t%s\n",$0,$2/$3,sp)}' > $SP1''.nonunq.count

#cat $SP1''.nonunq.count | grep -v 'Unknown' | awk '{sum += $4} END {print sum}' >> count.nonunq.sum
cat $SP1''.unq.count | grep -v 'Unknown' | awk '{sum += $4} END {print sum}' >> count.unq.sum

cd ..

done < ./samplenames_map_stats

#CHANGE FOLDER
cd count;

#cat *.nonunq.count > All_nonunq.count
#sed -i '1iTEfam\tcount\tmapped\tfrac\tspecies' All_nonunq.count
cat *.unq.count > All_unq.count
sed -i '1iTEfam\tcount\tmapped\tfrac\tspecies' All_unq.count

cd ..


printf "****DONE****\n"
