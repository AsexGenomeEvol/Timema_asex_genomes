#!/bin/bash
#
# is SP
# is a reference genome
# is the assembly folder

SP=$1
GENOME=$(basename $2 .fa)
ASM_DIR=$3

if [[ ! -s $2 ]]
then
	>&2 echo GENOME $2 IS MISSING;
	exit 1
fi

# bsub <<< """
bsub <<< """
#BSUB -L /bin/bash
#BSUB -J quast$SP
#BSUB -q bgee
#BSUB -o $3_quast.log
#BSUB -e $3_quast.err
#BSUB -n 32
#BSUB -M 64000000
#BSUB -R \"rusage[tmp=10000] span[ptile=32]\"

module add UHTS/Quality_control/quast/4.1

LOCAL_DIR=/scratch/local/monthly/kjaron/"$SP"_"$ASM_DIR"_quast
mkdir -p \$LOCAL_DIR

# copy required files
if [[ ! -s \$LOCAL_DIR/"$GENOME".fa ]]
then
	cp $2 \$LOCAL_DIR
fi

# go to local directory
cd \$LOCAL_DIR

# run quast
quast.py -e --no-sv --gene-finding -o "$GENOME"_quast -t 32 \
    -1 $TROOT/data/$SP/trimmed_reads/"$SP"_R1t_is_350.fq.gz  \
    -2 $TROOT/data/$SP/trimmed_reads/"$SP"_R2t_is_350.fq.gz \
    $GENOME.fa

mv "$GENOME"_quast $TROOT/data/$SP/assembly/$ASM_DIR/
rm "$GENOME".fa
rmdir \$LOCAL_DIR
"""
