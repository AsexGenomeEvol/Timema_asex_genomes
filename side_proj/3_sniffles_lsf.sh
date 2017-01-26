#!/bin/bash
#
# the first argument is the name of vcf output (and job)
# the second argument is a bam file with aligned long reads

bsub <<< """
#BSUB -L /bin/bash
#BSUB -J $1_snif
#BSUB -q dee
#BSUB -o $1.out
#BSUB -e $1.err
#BSUB -n 16
#BSUB -M 4000000
#BSUB -R \"rusage[tmp=30000] span[ptile=16]\"

WORKING_DIR=`pwd`

# make local directory for computations
mkdir -p /scratch/local/monthly/kjaron/$1
cd /scratch/local/monthly/kjaron/$1

# copy required files
cp $2 .

# run sniffles
sniffles -m `basename $2` -s 3 -t 16 -v $1.vcf

rm -f `basename $2`
mv * \$WORKING_DIR
rmdir /scratch/local/monthly/kjaron/$1
"""
