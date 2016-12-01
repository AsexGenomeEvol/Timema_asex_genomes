#!/bin/bash
#
# the first argument is the name of mapping output (and job)
# the second argument is an reference genome
# third  argument is file with reads to be mapped
# note that the output is rather big bam file, which will be compied to folder used for exection of the script

bsub <<< """
#BSUB -L /bin/bash
#BSUB -J $1_lr_map
#BSUB -q dee-hugemem
#BSUB -o $1.log
#BSUB -e $1.err
#BSUB -n 32
#BSUB -M 53554432
#BSUB -R \"rusage[tmp=70000] span[ptile=32]\"

module add UHTS/Analysis/samtools/1.3

WORKING_DIR=`pwd`

# make local directory for computations
mkdir -p /scratch/local/monthly/kjaron/$1
cd /scratch/local/monthly/kjaron/$1

# copy required files
cp $2 .
cp $3 .

# run bwa mem
ngmlr -t 31 -r `basename $2` -q `basename $3` | samtools view -bS - | samtools sort - > $1.bam

rm -f `basename $2` `basename $3` *.ngm
mv * \$WORKING_DIR
rmdir /scratch/local/monthly/kjaron/$1
"""
