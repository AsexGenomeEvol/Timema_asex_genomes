#!/bin/bash

# script that will generate apropriate script for abyss assembly
# and appropriate folder name for the assembly (printed on stdout)
# using template: $TROOT/C_contig_assembly/template_abyss_exec.sh

# 1. arg is path stick insect reads
# 2. is kmer
# then flags

USAGE="bash make_abyss_exec.sh SP KMER [--fewdata --fasteris --nose --nomse --nompe]"

if [ 2 -gt $# ]; then
  >&2 echo "invalid number of parameters, script usage:"
  >&2 echo $USAGE
  exit 1
fi

SP=$1
KMER=$2
# shift will remove $1 and $2 from variable containg all input parameters: $@
shift 2
# echo

DATAPATH=$TROOT/data/$SP/trimmed_reads
TARGETFOLDER="abyss_k$KMER"

TEMPNAME=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 16 | head -n 1)

sed "s/SP/$SP/g" $TROOT/C_contig_assembly/template_abyss_exec.sh > $TEMPNAME
sed -i "s/KMER/$KMER/g" $TEMPNAME

for i in "$@"
do
case $i in
## --fewdata : use also data from Dec 2016, (~15x of 350 pair-end)
    --fewdata)
      sed -i 's/pec //g' $TEMPNAME
      sed -i 's/pse_run2 //g' $TEMPNAME
      sed -i '/_run2/d' $TEMPNAME
      TARGETFOLDER+="_fewdata"
    ;;
## --fasteris : use mate pairs from fasteris (more conservative)
# maybe, just maybe, I should also modify the flag (FR vs RF stuff)
    --fasteris)
      # change path to fasteris MP instead of NxTrim
      sed -i 's/mp_nxtrim_FR/mp_fasteris_FR/g' $TEMPNAME
      # CHANGE OF ORIENTATION???
      TARGETFOLDER+="_fasteris"
    ;;
## --nose : no single end reads will be used
    --nose)
      # delete lines with pe and se library
      sed -i '/se/d' $TEMPNAME
      TARGETFOLDER+="_nose"
    ;;
## --nomse : only single end reads from pair end libraries will be used
    --nomse)
      sed -i 's/ mse//g' $TEMPNAME
      sed -i '/_se_mp/d' $TEMPNAME
      TARGETFOLDER+="_nomse"
    ;;
## --nompe : pair-end reads retrieved from mate-pair library won't be used
    --nompe)
      # delete lines starting by START_MPE ending by END_MPE which are flags around the library
      sed -i 's/pea //g' $TEMPNAME
      sed -i '/_is_225/d' $TEMPNAME
      TARGETFOLDER+="_nompe"
    ;;
## --filtered : pair-end reads mapping to contamination won't be used
    --filtered)
    # change file names to contaminaiton-filtered reads
      sed -i 's/fq.gz/no_contaminant.fastq.gz/g' $TEMPNAME
      TARGETFOLDER+="_filtered"
    ;;
    *)
      >&2 echo "Warning: unknown flag " ${i#*} " Usage:"
      >&2 echo $USAGE
    ;;
esac
done

mv $TEMPNAME $TARGETFOLDER.sh

echo $TARGETFOLDER
