#!/bin/bash
#
# the first argument is the name of mapping output (and job)
# the second argument is an indexed reference genome
# third and fourth arguments are reads to be mapped
# note that the output is rather bith bam file, which will be compied to folder used for exection of the script

bsub <<< """
#BSUB -L /bin/bash
#BSUB -J $1_map
#BSUB -q dee-hugemem
#BSUB -o $1.out
#BSUB -e $1.err
#BSUB -n 32
#BSUB -M 33554432
#BSUB -R \"rusage[tmp=70000] span[ptile=32]\"

module add UHTS/Aligner/bwa/0.7.13
module add UHTS/Analysis/samtools/1.3

WORKING_DIR=`pwd`

# make local directory for computations
mkdir -p /scratch/local/monthly/kjaron/$1
cd /scratch/local/monthly/kjaron/$1

# copy required files
cp $2* .
cp $3 .
cp $4 .

# run bwa mem
bwa mem -M -t 15 `basename $2` `basename $3` `basename $4` | samtools view -bS - | samtools sort - > $1.bam

mv $1.bam \$WORKING_DIR
rm -f `basename $2`* `basename $3` `basename $4`
rmdir /scratch/local/monthly/kjaron/$1
"""
