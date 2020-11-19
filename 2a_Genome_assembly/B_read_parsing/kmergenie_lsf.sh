#!/bin/bash
# 1. argument should be name of stick insect (5_Tge for instance)
# --filtered

FLAGS="-l 59 -k 89 -s 6 -t 16"

if [ "$2" = "--filtered" ] ; then
    SUFFIX=".no_contaminant.fastq.gz"
    DIR="kmergenie_filt"
    FLAGS="$FLAGS --diploid"
else
    SUFFIX="0.fq.gz"
    DIR="kmergenie"
fi

# uncomment
bsub <<< """
#BSUB -L /bin/bash
#BSUB -J $1_kmerg
#BSUB -q bgee
#BSUB -o $1_$DIR.out
#BSUB -e $1_$DIR.err
#BSUB -n 16
#BSUB -M 64000000
#BSUB -R \"rusage[tmp=80000]\"

# no need of kmergenie module, it is installed in /home/kjaron/bin/
module add R/latest

# make local directory for computations
mkdir -p /scratch/local/monthly/$USER/$1/$DIR
cd /scratch/local/monthly/$USER/$1/$DIR

cp \$TROOT/data/$1/trimmed_reads/*$SUFFIX .

ls *q.gz > read_files.list
kmergenie read_files.list $FLAGS

rm -f *$SUFFIX
mv /scratch/local/monthly/$USER/$1/$DIR \$TROOT/data/$1/trimmed_reads/

rmdir /scratch/local/monthly/$USER/$1
"""
