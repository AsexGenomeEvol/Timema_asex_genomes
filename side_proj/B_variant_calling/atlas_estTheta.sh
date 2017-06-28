#!/bin/bash
#
# is SP
# is a reference genome (VERSION)
# individual
# window_size

SP=$1
REF=$2
IND=$3
WINDOW=$4

BAMBASE="$IND"_to_"$REF"
BAM="$BAMBASE".bam
BASENAME="$IND"_to_"$REF"_w"$WINDOW"
LOCAL_DIR=/scratch/local/monthly/kjaron/"$SP"_"$BAM"_atlas_theta_w"$WINDOW"

if [[ ! -s $VROOT/data/$SP/mapping/$BAM ]]
then
	>&2 echo BAM FILE $BAM IS MISSING
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

# make local directory for computations
mkdir -p $LOCAL_DIR

# copy genome and mapped long insert size lib
cd $LOCAL_DIR
ln -s $VROOT/data/$SP/mapping/$BAM* .

atlas task=estimateTheta \
	bam=$BAM \
	window=$WINDOW \
	suppressWarnings verbose \
	1> "$BASENAME"_theta.log

mkdir -p $VROOT/data/$SP/variant_calls/$3/atlas
mv "$BAMBASE"_theta_estimates.txt \
	$VROOT/data/$SP/variant_calls/$IND/atlas/"$BASENAME"_theta_estimates.txt
mv "$BASENAME"_theta.log $VROOT/data/$SP/variant_calls/$IND/atlas/

rm $BAM*
rmdir \$LOCAL_DIR
"""
