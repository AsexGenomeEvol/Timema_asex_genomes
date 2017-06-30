#!/bin/bash
#
# is SP
# is a reference genome (VERSION)
# individual
# window_size

SP=$1
REF=$2
WINDOW=$3

BAMBASE=ref_to_"$REF"
BAM=$BAMBASE.bam
BAM1=ref_is350_to_"$REF".bam
BAM2=ref_is550_to_"$REF".bam
BAM3=ref_is700_to_"$REF".bam
BASENAME=ref_to_"$REF"_w"$WINDOW"
LOCAL_DIR=/scratch/local/monthly/kjaron/"$SP"_"$BAM"_atlas_theta_w"$WINDOW"

if [[ ! -s $TROOT/data/$SP/mapping/$BAM1 ]]
then
	>&2 echo BAM FILE $BAM1 IS MISSING
	exit 1
fi

bsub <<< """
#BSUB -L /bin/bash
#BSUB -J theta_$SP
#BSUB -q bgee
#BSUB -o job_"$BASENAME"_theta.log
#BSUB -e job_"$BASENAME"_theta.err
#BSUB -n 1
#BSUB -M 120000000

module add UHTS/Analysis/samtools/1.3

# make local directory for computations
mkdir -p $LOCAL_DIR

# copy genome and mapped long insert size lib
cd $LOCAL_DIR
ln -s $TROOT/data/$SP/mapping/$BAM1 .
ln -s $TROOT/data/$SP/mapping/$BAM2 .
ln -s $TROOT/data/$SP/mapping/$BAM3 .

samtools merge $BAM $BAM1 $BAM2 $BAM3
samtools index $BAM

atlas task=estimateTheta \
	bam=$BAM \
	window=$WINDOW \
	suppressWarnings verbose \
	1> "$BASENAME"_theta.log

mkdir -p $VROOT/data/$SP/variant_calls/ref/atlas
mv "$BAMBASE"_theta_estimates.txt \
	$VROOT/data/$SP/variant_calls/ref/atlas/"$BASENAME"_theta_estimates.txt
mv "$BASENAME"_theta.log $VROOT/data/$SP/variant_calls/ref/atlas/

rm $BAM*
rmdir \$LOCAL_DIR
"""
