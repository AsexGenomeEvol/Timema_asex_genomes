#!/bin/bash
#
# args:
# 1. name if stick insesct (3_Tms)
#
# input
# /scratch/beegfs/monthly/kjaron/timema_raw_reads/$1/is_[35]000/*R[12].fastq.gz
# output will be saved to
# /scratch/beegfs/monthly/kjaron/timema_trimmed/$1/mp_nxtrim_RF/$1_R[12]t_is_[35]000.fq.gz
# and
# /scratch/beegfs/monthly/kjaron/timema_trimmed/$1/$1_R[12]t_is_225.fq.gz and $1_se_mp.fq.gz
# ex: 1_Tdi_R1t_is_225

# testing if mate pairs are present

echo "Going to parse following files: " $TROOT/data/$1/raw_reads/is_[35]000/*R[12].fastq.gz

bsub <<< """
#BSUB -L /bin/bash
#BSUB -J $1_nxtrim
#BSUB -q bgee
#BSUB -o $1_mp_trim.out
#BSUB -e $1_mp_trim.err
#BSUB -n 4
#BSUB -M 4194304
#BSUB -R \"rusage[tmp=40000] span[ptile=4]\"

READSPATH=$TROOT/data/$1/raw_reads
TARGETPATH=$TROOT/data/$1/trimmed_reads

# make local directory for computations
mkdir -p /scratch/local/monthly/$USER/$1
cd /scratch/local/monthly/$USER/$1

# copy required files
cp \$READSPATH/is_[35]000/*R[12].fastq.gz .

# is_3000 libs have L004 in name, is_5000 have L005
nxtrim -1 *L004*R1* -2 *L004*R2* -O is_3000 --preserve-mp --separate &
nxtrim -1 *L005*R1* -2 *L005*R2* -O is_5000 --preserve-mp --separate &

wait

zcat is_3000_R1.mp.fastq.gz is_3000_R1.unknown.fastq.gz | gzip > $1_R1t_is_3000.fq.gz &
zcat is_3000_R2.mp.fastq.gz is_3000_R2.unknown.fastq.gz | gzip > $1_R2t_is_3000.fq.gz &
zcat is_5000_R1.mp.fastq.gz is_5000_R1.unknown.fastq.gz | gzip > $1_R1t_is_5000.fq.gz &
zcat is_5000_R2.mp.fastq.gz is_5000_R2.unknown.fastq.gz | gzip > $1_R2t_is_5000.fq.gz &

wait

cat is_3000_R1.pe.fastq.gz is_5000_R1.pe.fastq.gz > PE_R1t.fastq.gz &
cat is_3000_R2.pe.fastq.gz is_5000_R2.pe.fastq.gz > PE_R2t.fastq.gz &

wait

# delete input files (they are no needed anymore)
rm -f *.fastq.gz &

gzip *_se_mp.fq &
gzip *_R1t_is_225.fq &
gzip *_R2t_is_225.fq &

wait

# MOVE BACK RESULTS
mkdir -p \$TARGETPATH/mp_nxtrim_FR/
mv *_is_[35]000.fq.gz \$TARGETPATH/mp_nxtrim_FR/
mv *_se_mp.fq.gz *_R1t_is_225.fq.gz *_R2t_is_225.fq.gz \$TARGETPATH/

#rmdir /scratch/local/monthly/$USER/$1
"""
