#!/bin/bash

clear
while read line; do
  read -a twoSp <<< $line
  printf "\n****processing ${twoSp[0]}*****\n"
  SP1=${twoSp[0]}

printf "**running $SP1**\n\n"

source /scratch/beegfs/monthly/ptranvan/Software/PASTEClassifier/1.0/config/setEnv.sh
module add SequenceAnalysis/Repeat/trf/4.07b;
module add SequenceAnalysis/HMM-Profile/hmmer/3.1b2;


mkdir $SP1''_class;

cd $SP1''_class;

cat ../PASTEClassifier_parallelized.cfg | sed s"/project_name: project_name/project_name: $SP1\_class/"g > tmp
cat tmp | sed s"/project_dir: \/project_dir/project_dir: \/scratch\/beegfs\/monthly\/jbast\/TEclassi\/$SP1\_class/"g > PASTEClassifier_parallelized.cfg
rm tmp

ln -s /scratch/beegfs/monthly/ptranvan/Software/REPET/repbase/RepBase20.05_REPET.embl/repbase20.05_ntSeq_cleaned_TE.fa repbase20.05_ntSeq_cleaned_TE.fa
ln -s /scratch/beegfs/monthly/ptranvan/Software/REPET/repbase/RepBase20.05_REPET.embl/repbase20.05_aaSeq_cleaned_TE.fa repbase20.05_aaSeq_cleaned_TE.fa
ln -s /scratch/beegfs/monthly/ptranvan/Software/REPET/repbase/ProfilesBankForREPET_Pfam27.0_GypsyDB.hmm ProfilesBankForREPET_Pfam27.0_GypsyDB.hmm
ln -s /scratch/beegfs/monthly/jbast/TEclassi/HG/$SP1''_b3v06.max.transcripts.func.fasta HG.fasta
ln -s /scratch/beegfs/monthly/jbast/TEclassi/rDNA/rDNA_silva.fa rDNA_silva.fa

cp /scratch/beegfs/monthly/jbast/masking_merged_Timema/raw_libraries/$SP1''_merged.TElib.fa.centroid95min500.cvt ./

cd ..

done < samplenames


printf "****DONE****\n"
