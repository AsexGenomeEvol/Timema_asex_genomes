#!/bin/bash

#get variables to use in bash from two columns of a file per row
#map reads to timema genomes

#folder needed:
# reads
# map
# genomes


clear
while read line; do
  read -a twoSp <<< $line
  printf "\n****processing ${twoSp[0]}*****\n"
  SP1=${twoSp[0]}

#get reads
printf "**get and process reads for $SP1**\n\n"
#cp /scratch/beegfs/monthly/kjaron/timema_assembly/data/$SP1''/trimmed_reads/is_350/$SP1''_is_350_SND393_L001_R1t.fq.gz ./reads
#cp /scratch/beegfs/monthly/kjaron/timema_assembly/data/$SP1''/trimmed_reads/is_350/$SP1''_is_350_SND393_L001_R2t.fq.gz ./reads
#cd reads/
#gunzip *.gz
#cd ..

#rename read headers
cd reads
#sed 's/ 1:N:0:\(.*\)$/\#\1\/1/g' $SP1''_is_350_SND393_L001_R1t.fq > $SP1''_is_350_SND393_L001_R1t_rn.fq
#sed 's/ 2:N:0:\(.*\)$/\#\1\/2/g' $SP1''_is_350_SND393_L001_R2t.fq > $SP1''_is_350_SND393_L001_R2t_rn.fq
cd ..

#index genomes
cd genomes;
#for f in *.fasta; do bwa index $f; done
cd ..

#mapping (-t is number of CPUs)
printf "**mapping $SP1**\n\n"
bwa mem -t 60 ./genomes/$SP1''_b3v08.fasta ./reads/$SP1''_is_350_SND393_L001_R1t_rn.fq ./reads/$SP1''_is_350_SND393_L001_R2t_rn.fq > ./map/$SP1''_aln.sam;

#samtools
printf "**samtools $SP1**\n\n"
cd map
samtools view -@ 60 -bS -F 4 $SP1''_aln.sam -o $SP1''_aln.bam
samtools sort -@ 40 $SP1''_aln.bam -o $SP1''_aln.sorted.bam
samtools index $SP1''_aln.sorted.bam
#samtools flagstat $SP1''_aln.sorted.bam &> $SP1''_aln.sorted.bam.stat; done
cd ..

#remove unnecessary files
  rm ./map/$SP1''_aln.sam;
  rm ./map/$SP1''_aln.bam;

#HTSeq-count
printf "**HTseq-count $SP1**\n\n"
htseq-count -f bam -r name -s no -t similarity -i Target --nonunique none map/$SP1''_aln.sorted.bam gff/$SP1''_b3v08.fa.out.len80.gff &> ./count/$SP1''_htseq_unq.cnt

done < ./samplenames

#module rm UHTS/Aligner/bwa/0.7.13;
#module rm UHTS/Analysis/samtools/1.3;
#module rm Development/java/1.8.0_202;
#module rm UHTS/Analysis/HTSeq/0.9.1;


printf "****DONE****\n"
