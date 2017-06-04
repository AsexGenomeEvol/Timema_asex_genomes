#!/bin/bash
#
# is SP
# is a reference genome (VERSION)
# individual
# window_size

SP=$1
BAM="$3"_to_"$2".bam
WINDOW=$4

if [[ ! -s $TROOT/data/$SP/mapping/$BAM ]]
then
	>&2 echo BAM FILE $BAM IS MISSING
	exit 1
fi

bsub <<< """
#BSUB -L /bin/bash
#BSUB -J theta_$1
#BSUB -q bgee
#BSUB -o job_"$3"_to_"$2"_w"$WINDOW"_theta.log
#BSUB -e job_"$3"_to_"$2"_w"$WINDOW"_theta.err
#BSUB -n 1
#BSUB -M 20000000

# make local directory for computations
LOCAL_DIR=/scratch/local/monthly/kjaron/"$SP"_"$BAM"_atlas_theta
mkdir -p \$LOCAL_DIR

# copy genome and mapped long insert size lib
cd \$LOCAL_DIR
ln -s $TROOT/data/$SP/mapping/$BAM* .

python3 $TROOT/N_variant_calling/create_BED_withmin_window_size.py \
	$SP $3 $WINDOW > relevant_windows.bed
atlas task=estimateTheta bam=$BAM suppressWarnings window=relevant_windows.bed \
	1> "$3"_to_"$2"_w"$WINDOW"_theta.log 2> "$3"_to_"$2"_w"$WINDOW"_theta.err

mkdir -p $TROOT/data/$SP/variant_calls/$3/atlas
mv $(basename $BAM .bam)_theta_estimates.txt \
	$TROOT/data/$SP/variant_calls/$3/atlas/$(basename $BAM .bam)_w"$WINDOW"_theta_estimates.txt
mv "$3"_to_"$2"_w"$WINDOW"_theta.* $TROOT/data/$SP/variant_calls/$3/atlas/

rm $BAM*
rmdir \$LOCAL_DIR
"""
