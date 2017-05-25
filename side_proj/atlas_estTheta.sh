#!/bin/bash
#
# is SP
# is a reference genome (VERSION)
# individual

SP=$1
BAM="$3"_to_"$2".bam

if [[ ! -s $TROOT/data/$SP/mapping/$BAM ]]
then
	>&2 echo BAM FILE $BAM IS MISSING
	exit 1
fi

bsub <<< """
#BSUB -L /bin/bash
#BSUB -J theta_$1
#BSUB -q bgee
#BSUB -o "$3"_to_"$2"_theta.out
#BSUB -e "$3"_to_"$2"_theta.err
#BSUB -n 1
#BSUB -M 20000000

# make local directory for computations
LOCAL_DIR=/scratch/local/monthly/kjaron/"$SP"_"$BAM"_atlas_theta
mkdir -p \$LOCAL_DIR

# copy genome and mapped long insert size lib
cd \$LOCAL_DIR
ln -s $TROOT/data/$SP/mapping/$BAM* .

atlas task=estimateTheta bam=$BAM suppressWarnings window=1000

mkdir -p $TROOT/data/$SP/variant_calls/$3/atlas
mv $(basename $BAM .bam)_theta_estimates.txt $TROOT/data/$SP/variant_calls/$3/atlas/

rm $BAM*
rmdir \$LOCAL_DIR
"""
