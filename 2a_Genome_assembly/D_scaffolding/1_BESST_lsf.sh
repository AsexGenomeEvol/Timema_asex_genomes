#!/bin/bash
#
# the first argument is the SP
# the second argument is the folder assembly
# third argument is name of reference

SP=$1
GENOME=$3

OUTPUT=$2_BESST

LOCAL_DIR=/scratch/local/monthly/kjaron/$1$2
ASM_DIR=\$TROOT/data/$SP/assembly/$2

bsub <<< """
#BSUB -L /bin/bash
#BSUB -J Bscf_$1
#BSUB -q bgee
#BSUB -o $OUTPUT.log
#BSUB -e $OUTPUT.err
#BSUB -n 1
#BSUB -M 33554432
#BSUB -R \"rusage[tmp=50000] span[ptile=1]\"

# make local directory for computations
mkdir -p $LOCAL_DIR/temp
export TMPDIR=$LOCAL_DIR/temp
cd $LOCAL_DIR

ln -s $ASM_DIR/$GENOME .
ln -s $ASM_DIR/mapping/*bam .
ln -s $ASM_DIR/mapping/*bam.bai .

# run BESST script
python ~/src/BESST/runBESST -c $GENOME \
  -f $SP""_350.bam $SP""_550.bam $SP""_700.bam $SP""_3000.bam $SP""_5000.bam \
  -o $OUTPUT -orientation fr fr fr fr fr

mv $OUTPUT $TROOT/data/$SP/assembly/
"""
