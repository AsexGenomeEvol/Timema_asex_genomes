#!/bin/bash
# 1. argument should be name of stick insect (5_Tge for instance)
# 2. argument should be kmer for assembly

USAGE="bash 1_abyss.sh SP KMER [--fewdata --fasteris --nose --nomse --nompe --filtered]"

if [ 2 -gt $# ]; then
  >&2 echo "invalid number of parameters, script usage:"
  >&2 echo $USAGE
  exit 1
fi

SP=$1
KMER=$2
shift 2


# $TROOT/B_read_parsing/verify_trimming.sh $SP

if [[ $? -ne 0 ]]
then
	exit 1
fi

ASM_DIR=$($TROOT/C_contig_assembly/make_abyss_exec.sh $SP $KMER $@)

if [[ -s $TROOT/data/$SP/assembly/$ASM_DIR ]]
then
	echo ERROR: assembly $ASM_DIR was already performed
	echo check $TROOT/data/$SP/assembly/$ASM_DIR
	rm $ASM_DIR.config
	exit 1
fi

echo "tests passed. Submiting job..."

# uncomment
bsub <<< """
#BSUB -L /bin/bash
#BSUB -J $SP""abys$KMER
#BSUB -q bgee
#BSUB -o $ASM_DIR.out
#BSUB -e $ASM_DIR.err
#BSUB -n 32
#BSUB -M 154217728
#BSUB -R \"rusage[tmp=30000] span[ptile=32]\"

# no need of module, abyss is installed in /home/kjaron/bin/

LOCALDIR=/scratch/local/monthly/kjaron/$SP/$ASM_DIR

# make local directory for computations
mkdir -p \$LOCALDIR/temp
export TMPDIR=\$LOCALDIR/temp

mv $ASM_DIR.sh \$LOCALDIR

cd \$LOCALDIR
ln -s $TROOT/data/$SP/trimmed_reads/* .

bash $ASM_DIR.sh

mkdir -p $TROOT/data/$SP/assembly/
mv \$LOCALDIR $TROOT/data/$SP/assembly/
"""
