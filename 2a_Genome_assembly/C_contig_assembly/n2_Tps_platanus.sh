#!/bin/bash
#
# script description
# text in <> should be replaced according the job
# make one folder for every job you sumbit, copy, move or link there input files
# modify and submit this script

#BSUB -L /bin/bash
#BSUB -J Tps_planatus
#BSUB -q dee-hugemem
#BSUB -o Tps_platanus.out
#BSUB -e Tps_platanus.err
#BSUB -n 32
#BSUB -M 256000
#BSUB -R "rusage[tmp=10000] span[ptile=32]"

module add UHTS/Assembler/platanus/1.2.1

WORKING_DIR=`pwd`
# make local directory for computations
mkdir -p /scratch/local/monthly/plat
cd /scratch/local/monthly/plat

# copy all input files from working directory to local disc
cp $WORKING_DIR/Tps_R[12]t_is_[2357][250][50].fq .
INPUT_FILES=`ls`

Platanus assemble -o Tps_plt -f ./Tps_R[12]t_is_[2357][250][50].fq -t 32 -m 256

# remove input files
rm -f $INPUT_FILES
# move results to your scratch
mv * $WORKING_DIR
# remove directory
rmdir /scratch/local/monthly/plat
