#!/bin/bash

#get variables to use in bash from two columns of a file per row

#get reseq trimmed reads from cp /archive/dee/schwander/kjaron/timema_reseq_reads.tar ./


#map reads to timema genomes v7

#folder needed:
# reads
# map
# genomes


########
#BEFORE#
########

#before, make filelist
#first remove all se reads (grep -v), then:
#cat filelist | grep -v '_se.' | sed 's/\/scratch\/beegfs\/monthly\/odegaspe\/timema_reseq_trim\///g' | sed 's/_R1.cleaned.fastq.gz//g' | sed 's/_R2.cleaned.fastq.gz//g' | sed 's/_OBIWAN/\t_OBIWAN/g' | sed 's/1_Tdi/1_Tdi\t/g' | sed 's/1_Tps/1_Tps\t/g' | sed 's/2_Tcm/2_Tcm\t/g' | sed 's/2_Tsi/2_Tsi\t/g' | sed 's/3_Tce/3_Tce\t/g' | sed 's/3_Tms/3_Tms\t/g' | sed 's/4_Tbi/4_Tbi\t/g' | sed 's/4_Tte/4_Tte\t/g' | sed 's/5_Tge/5_Tge\t/g' | sed 's/5_Tpa/5_Tpa\t/g' sed -e 's/_L/\t_L/g' > filelist.txt 

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
#  SP3=${twoSp[2]}
#  SP4=${twoSp[3]}


printf "\n**processing $SP1$SP2**\n"

#get reads
#cp /scratch/beegfs/monthly/odegaspe/timema_reseq_trim/$SP1$SP2$SP3$SP4''_R1.cleaned.fastq.gz ./reads
#cp /scratch/beegfs/monthly/odegaspe/timema_reseq_trim/$SP1$SP2$SP3$SP4''_R2.cleaned.fastq.gz ./reads
#cd reads/
#gunzip *.gz
#cd ..

#rename read headers
#cd reads
#zcat ../reads_reseq_trimmed/$SP1$SP2$SP3$SP4''_R1.cleaned.fastq.gz | sed 's/ 1:N:0:\(.*\)$/\#\1\/1/g' > $SP1$SP2$SP3$SP4''_R1.cleaned_rn.fastq
#zcat ../reads_reseq_trimmed/$SP1$SP2$SP3$SP4''_R2.cleaned.fastq.gz | sed 's/ 2:N:0:\(.*\)$/\#\1\/2/g' > $SP1$SP2$SP3$SP4''_R2.cleaned_rn.fastq
#cd ..


#merge reads
#cd reads
#cat $SP1$SP2''*''_R1.cleaned_rn.fastq > $SP1$SP2''_R1.cleaned_rn_mrg.fastq
#cat $SP1$SP2''*''_R2.cleaned_rn.fastq > $SP1$SP2''_R2.cleaned_rn_mrg.fastq
#cd ..


#index genomes
#cd genomes;
#for f in *.fa; do bwa index $f; done
#cd ..

#mapping (-t is number of CPUs)
printf "**mapping $SP1**\n\n"
bwa mem -t 60 ./genomes/$SP1''_b3v08.fasta ./reads/$SP1$SP2''_R1.cleaned_rn_mrg.fastq.gz ./reads/$SP1$SP2''_R2.cleaned_rn_mrg.fastq.gz > ./map/$SP1$SP2''_aln.sam;

#samtools

printf "**samtools $SP1**\n\n"
cd map
samtools view -@ 60 -bS -F4 $SP1$SP2''_aln.sam -o $SP1$SP2''_aln.bam
rm ./$SP1$SP2''_aln.sam
samtools sort -@ 60 $SP1$SP2''_aln.bam -o $SP1$SP2''_aln.sorted.bam
samtools index $SP1$SP2''_aln.sorted.bam
#samtools flagstat $SP1$SP2''_aln.sorted.bam &> $SP1$SP2''_aln.sorted.bam.stat; done
cd ..

#remove unnecessary files
  #rm ./map/$SP1$SP2''_aln.sam;
  rm ./map/$SP1$SP2''_aln.bam;

#HTSeq-count
printf "**HTseq-count $SP1**\n\n"
#htseq-count -f bam -r name -s no -t similarity -i Target --nonunique none map/$SP1$SP2''_aln.sorted.bam gff/$SP1''_b3v08.fa.out.len80.gff &> ./count/$SP1$SP2''_htseq_unq.cnt


done < ./filelist5.txt


printf "****DONE****\n"
