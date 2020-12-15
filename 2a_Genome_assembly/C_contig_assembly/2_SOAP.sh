#!/bin/bash
# 1. argument should be name of stick insect (5_Tge for instance)
# 2. argument should be kmer for assembly

USAGE="bash 2_SOAP.sh SP KMER [--fewdata --fasteris --nose --nomse --nompe]"

if [ 2 -gt $# ]; then
  >&2 echo "invalid number of parameters, script usage:"
  >&2 echo $USAGE
  exit 1
fi

SP=$1
KMER=$2
shift 2

$TROOT/B_read_parsing/verify_trimming.sh $SP

if [[ $? -ne 0 ]]
then
	exit 1
fi

ASM_DIR=$($TROOT/C_contig_assembly/make_SOAP_setting.sh $SP $KMER $@)

if [[ -s $TROOT/data/$SP/assembly/$ASM_DIR ]]
then
	echo ERROR: assembly $ASM_DIR was already performed
	echo check $TROOT/data/$SP/assembly/$ASM_DIR
	rm $ASM_DIR.config
	exit 1
fi

echo "tests passed. Submiting job..."
echo "assembly directory:" $ASM_DIR

# uncomment
bsub <<< """
#BSUB -L /bin/bash
#BSUB -J $SP""SOAP$KMER
#BSUB -q bgee
#BSUB -o $ASM_DIR.log
#BSUB -e $ASM_DIR.err
#BSUB -n 32
#BSUB -M 157286400
#BSUB -R \"rusage[tmp=40000] span[ptile=32]\"
# variable COMPUTER can contain string specifing computer to compute on
$COMPUTER

# no need of module, SOAP is installed in /home/kjaron/bin/
LOCALDIR=/scratch/local/monthly/kjaron/$SP/$ASM_DIR

# make local directory for computations
mkdir -p \$LOCALDIR/temp
export TMPDIR=\$LOCALDIR/temp

mv $ASM_DIR.config \$LOCALDIR
cd \$LOCALDIR

OUT_BASE=$SP""_SOAP_k$KMER
# if K > 63 it should run the second exectuble
SOAPdenovo-63mer all -s $ASM_DIR.config -p 32 -a 120 -K $KMER -R -o \$OUT_BASE

# # compute directly assembly stats (does not work for reason I do not get)
# $TROOT/scripts/write_stats.sh \$OUT_BASE"".contig > \$OUT_BASE""_ctg.stats &
# $TROOT/scripts/write_stats.sh \$OUT_BASE"".scafSeq > \$OUT_BASE""_scf.stats

wait

rm -r temp
mkdir -p $TROOT/data/$SP/assembly/
mv \$LOCALDIR $TROOT/data/$SP/assembly/
"""
