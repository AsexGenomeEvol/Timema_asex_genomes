#!/bin/bash
#
# is SP
# is a reference genome (VERSION)
# individual

SP=$1
GENOME="$1"_"$2".fa
BAM="$3"_to_"$2".bam

if [[ ! -s $TROOT/data/$SP/reference/$GENOME ]]
then
	>&2 echo GENOME $2 IS MISSING IN THE REFERENCE FOLDER;
	exit 1
fi

if [[ ! -s $TROOT/data/$SP/mapping/$BAM ]]
then
	>&2 echo BAM FILE $BAM IS MISSING
	exit 1
fi

bsub <<< """
#BSUB -L /bin/bash
#BSUB -J callMLE_$1
#BSUB -q bgee
#BSUB -o "$3"_to_"$2"_callMLE.out
#BSUB -e "$3"_to_"$2"_callMLE.err
#BSUB -n 1
#BSUB -M 20000000

module add UHTS/Analysis/samtools/1.3

# make local directory for computations
LOCAL_DIR=/scratch/local/monthly/kjaron/"$SP"_"$BAM"_atlas_variants
mkdir -p \$LOCAL_DIR

# copy genome and mapped long insert size lib
cd \$LOCAL_DIR
ln -s $TROOT/data/$SP/mapping/$BAM* .
ln -s $TROOT/data/$SP/reference/$GENOME* .

atlas task=callMLE bam=$BAM fasta=$GENOME vcf suppressWarnings

mkdir -p $TROOT/data/$SP/variant_calls/$3/atlas
mv $(basename $BAM .bam)_MLEGenotypes.vcf.gz $TROOT/data/$SP/variant_calls/$3/atlas/

rm $BAM*
rm $GENOME*

rmdir \$LOCAL_DIR
"""
