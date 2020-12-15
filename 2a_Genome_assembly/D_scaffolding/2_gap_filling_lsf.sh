#!/bin/bash
#
# the first argument is the name of stick insect (5_Tge for instance)
# he reference genome
# an assembly folder
# (optional) --filtered (to use contamination filtered reads)

SOAP_CONFIG_FLAGS="--fewdata --nomse --nompe"

if [ "$4" = "--filtered" ] ; then
    FILT_SUFFIX="_filtered"
    SOAP_CONFIG_FLAGS="$SOAP_CONFIG_FLAGS --filtered"
fi

SP=$1

TARGET=$(basename ${2%.*})_GC.fa
TARGET_DIR=$3

JOB_ID=$TARGET_DIR""_GC
LOCAL_DIR=/scratch/local/monthly/kjaron/$JOB_ID

bsub <<< """
#BSUB -L /bin/bash
#BSUB -J $JOB_ID
#BSUB -q bgee
#BSUB -o $JOB_ID.log
#BSUB -e $JOB_ID.err
#BSUB -n 64
#BSUB -M 500000000
#BSUB -R \"rusage[tmp=20000] span[ptile=64]\"

# make local directory for computations
mkdir -p $LOCAL_DIR/temp
export TMPDIR=$LOCAL_DIR/temp

cp $2 $LOCAL_DIR
cd $LOCAL_DIR

# script make_SOAP_setting.sh creates a setting file for SOAP tools
# setting file is used only to specify libraries for gap filling
CONFIG=$TROOT/C_contig_assembly/make_SOAP_setting.sh $SP 0 $SOAP_CONFIG_FLAGS $4

GapCloser \
  -a $(basename $2) -b $CONFIG.config \
  -o $TARGET -t 64 -l 150

mv $TARGET* $TROOT/data/$SP/assembly/$TARGET_DIR
rm $(basename $2)
"""
