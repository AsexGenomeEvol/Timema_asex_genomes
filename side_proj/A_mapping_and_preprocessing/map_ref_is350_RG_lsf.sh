#!/bin/bash
#
# the first argument is the name of species (1_Tdi)
# the second argument is a version indexed reference genome (b3v04)
# the third argument is read group
# (ex: "@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1" )
# note that the output is rather bith bam file, which will be compied to folder data/SP/mapping

BAM=ref_is350_to_$2.bam

if [[ -s $TROOT/data/$1/mapping/$BAM ]]
then
	>&2 echo MAPPING FILE $BAM EXISTS ALREADY;
	exit 1
fi

bsub <<< """
#BSUB -L /bin/bash
#BSUB -J "$1"_"$2"_map
#BSUB -q bgee
#BSUB -o ref_is350_to_"$2".log
#BSUB -e ref_is350_to_"$2".err
#BSUB -n 16
#BSUB -M 41943040
#BSUB -R \"rusage[tmp=30000] span[ptile=16]\"

module add UHTS/Aligner/bwa/0.7.13
module add UHTS/Analysis/samtools/1.3

# make local directory for computations
LOCAL_DIR=/scratch/local/monthly/kjaron/map_$1
mkdir -p \$LOCAL_DIR/temp
export TMPDIR=\$LOCAL_DIR/temp

cd \$LOCAL_DIR

# copy reference files
cp $TROOT/data/$1/reference/$1_$2* .


# run bwa mem
bwa mem -M -t 15 -R \"$3\" \
  $1_$2.fa \
  $TROOT/data/$1/trimmed_reads/$1_R1t_is_350.fq.gz \
  $TROOT/data/$1/trimmed_reads/$1_R2t_is_350.fq.gz \
  | samtools sort -O bam - > $BAM

mkdir -p $TROOT/data/$1/mapping
mv $BAM $TROOT/data/$1/mapping/

rm -rf temp
rm -f $1_$2*
rmdir \$LOCAL_DIR
"""
