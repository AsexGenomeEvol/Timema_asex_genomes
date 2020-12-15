#!/bin/bash
#
# the 1. argument is the name of mapping output (and job)
# the 2.  argument is a reference genome
# 3. and 4. arguments are reads to be mapped
# note that the output is rather bith bam file, which will be compied to folder used for exection of the script


# test for reference
if [[ ! -s $2 ]]
then
	echo REFERENCE $2 IS MISSING
	exit 1
fi

# test for reads R1
if [[ ! -s $3 ]]
then
	echo READ FILE $3 IS MISSING
	exit 1
fi

# test for reads R2
if [[ ! -s $4 ]]
then
	echo READ FILE $4 IS MISSING
	exit 1
fi

bsub <<< """
#BSUB -L /bin/bash
#BSUB -J map$1
#BSUB -q bgee
#BSUB -o $1.out
#BSUB -e $1.err
#BSUB -n 16
#BSUB -M 52428800
#BSUB -R \"rusage[tmp=30000] span[ptile=16]\"

module add UHTS/Aligner/bwa/0.7.15;
module add UHTS/Analysis/samtools/1.3

WORKING_DIR=`pwd`
LOCAL_DIR=/scratch/local/monthly/kjaron/$1

# make local directory for computations
mkdir -p \$LOCAL_DIR/temp
export TMPDIR=\$LOCAL_DIR/temp

# copy required files
cp $2* \$LOCAL_DIR &
cp $3 \$LOCAL_DIR &
cp $4 \$LOCAL_DIR

wait

cd \$LOCAL_DIR

# run bwa mem using options from BESST script
bwa mem -t 15 -w 0 -O 99 `basename $2` `basename $3` `basename $4` | samtools view -bS - | samtools sort - > $1.bam

mv $1* \$WORKING_DIR
rm `basename $2`* `basename $3` `basename $4`
"""
