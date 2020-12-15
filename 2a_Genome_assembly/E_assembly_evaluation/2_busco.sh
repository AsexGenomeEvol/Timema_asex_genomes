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

bsub <<< """
#BSUB -L /bin/bash
#BSUB -J busco$SP
#BSUB -q bgee
#BSUB -o $3_busco.log
#BSUB -e $3_busco.err
#BSUB -n 32
#BSUB -M 32000000
#BSUB -R \"rusage[tmp=10000] span[ptile=32]\"

export PATH=\$PATH:/scratch/beegfs/monthly/ptranvan/Software/busco/2.0/
module add Blast/ncbi-blast/2.2.31+
module add SequenceAnalysis/HMM-Profile/hmmer/3.1b2
module add SequenceAnalysis/GenePrediction/augustus/3.2.2
export AUGUSTUS_CONFIG_PATH=/scratch/beegfs/monthly/kjaron/augustus_config

LOCAL_DIR=/scratch/local/monthly/kjaron/"$SP"_"$ASM_DIR"_busco
mkdir -p \$LOCAL_DIR

# copy required files
cp $2 \$LOCAL_DIR

# go to local directory
cd \$LOCAL_DIR

# run busco
BUSCO.py -i $GENOME.fa -o "$GENOME"_busco -m geno \
         -l /scratch/beegfs/monthly/kjaron/busco_ref/insecta_odb9 -c 32


cp run_"$GENOME"_busco/short_summary* $TROOT/data/$SP/assembly/$ASM_DIR/
cp run_"$GENOME"_busco/full_table* $TROOT/data/$SP/assembly/$ASM_DIR/
cp run_"$GENOME"_busco/missing_busco* $TROOT/data/$SP/assembly/$ASM_DIR/

tar czf "$SP"_"$GENOME"_busco.tar.gz run_"$GENOME"_busco
mv "$SP"_"$GENOME"_busco.tar.gz $TROOT/data/$SP/assembly/$ASM_DIR/
rm -r run_"$GENOME"_busco
rm -r tmp
rm "$GENOME".fa
rmdir \$LOCAL_DIR
"""
