#!/bin/bash
#
# is SP
# is a reference genome (VERSION)
# individual

SP=$1
GENOME="$1"_"$2".fa
BAM="$3"_to_"$2".bam
DELLY="$3"_to_"$2".bcf

if [[ ! -s $TROOT/data/$SP/reference/$GENOME.fai ]]
then
	>&2 echo INDEX OF GENOME $2 IS MISSING IN THE REFERENCE FOLDER;
	exit 1
fi

if [[ ! -s $TROOT/data/$SP/mapping/$BAM.bai ]]
then
	>&2 echo INDEX OF BAM FILE $BAM IS MISSING
	exit 1
fi

#BSUB -L /bin/bash
#BSUB -J "$3"_delly
#BSUB -q bgee
#BSUB -o "$3"_delly.log
#BSUB -e "$3"_delly.err
#BSUB -n 1
#BSUB -M 200000000
#BSUB -R \"rusage[tmp=50000]\"

LOCAL_DIR=/scratch/local/monthly/kjaron/"$SP"_"$BAM"_delly
mkdir -p \$LOCAL_DIR/temp
export TMPDIR=\$LOCAL_DIR/temp

# copy genome and mapped long insert size lib
cp $TROOT/data/$SP/reference/$GENOME \
   $TROOT/data/$SP/reference/$GENOME.fai \$LOCAL_DIR
cp $TROOT/data/$SP/mapping/$BAM \
   $TROOT/data/$SP/mapping/$BAM.bai \$LOCAL_DIR

# go to local directory
cd \$LOCAL_DIR

delly call -g $GENOME $BAM -o $DELLY

mkdir -p $TROOT/data/$SP/variant_calls/$3/delly
mv $DELLY $TROOT/data/$SP/variant_calls/$3/delly/

rm $GENOME* $BAM*
rm -r temp

rmdir \$LOCAL_DIR
