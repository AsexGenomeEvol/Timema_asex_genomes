#!/bin/bash
#
# 1. argument is the folder of the species (ex. 1_Tdi or 2_Tge)
# 2. argument: ilbrary (ex. is_350)
# 3. (optional) is a PATTERN read files have to start, for instance for sequences from Fasteris from December 2016: (ex. 1612)
# 4. (optional) is a suffix added to output files

# trimmed files are moved to folder from where the script was executed
# reads are expected to be found /scratch/beegfs/monthly/kjaron/data/$1/raw_reads/$2
# file ~/adapters/AllIllumina-PEadapters.fa will be used for trimming adapters

# test for input
if [[ ! "$(ls -A $TROOT/data/$1/raw_reads/$2/)" ]]
then
	echo NO FILES TO TRIM AT: $TROOT/data/$1/raw_reads/$2/;
	exit 1
fi

# test for output
if [[ -s $TROOT/data/$1/trimmed_reads/$1""_R1t_$2$4.fq.gz ]]
then
	echo FILE $TROOT/data/$1/trimmed_reads/$1""_R1t_is_350.fq.gz ALREADY EXISTS;
	exit 1
fi

bsub <<< """
#BSUB -L /bin/bash
#BSUB -J $1_trmm
#BSUB -q bgee
#BSUB -o $1$2$3$4.log
#BSUB -e $1$2$3$4.err
#BSUB -n 16
#BSUB -M 25165824
#BSUB -R \"rusage[tmp=70000] span[ptile=16]\"

module add UHTS/Analysis/trimmomatic/0.36;

READSPATH=$TROOT/data/$1/raw_reads/$2
LOCALDIR=/scratch/local/daily/$USER/$1$2$3$4
export TMPDIR=\$LOCALDIR/temp
TARGETPATH=$TROOT/data/$1/trimmed_reads

mkdir -p \$LOCALDIR/temp
mkdir -p \$TARGETPATH

cd \$LOCALDIR
cp ~/adapters/AllIllumina-PEadapters.fa .

# concatinate all reads of same type, trim them and remove the concatined file (the source remains unchanged)
fR1=$1""_temp_R1.fq.gz
fR2=$1""_temp_R2.fq.gz
for read_file in \$(ls \$READSPATH/$3*R1.f*); do
    	cat \$read_file >> \$fR1
done
for read_file in \$(ls \$READSPATH/$3*R2.f*); do
       	cat \$read_file >> \$fR2
done

trimmomatic PE -threads 16 \$fR1 \$fR2 $1""_R1t_$2$4.fq.gz $1""_R1np_$2$4.fq.gz $1""_R2t_$2$4.fq.gz $1""_R2np_$2$4.fq.gz ILLUMINACLIP:AllIllumina-PEadapters.fa:3:25:6 LEADING:9 TRAILING:9 SLIDINGWINDOW:4:15 MINLEN:96

mv $1""_R[12]t_$2$4.fq.gz \$TARGETPATH
mv $1""_R[12]np_$2$4.fq.gz \$TARGETPATH

rm \$fR1
rm \$fR2

# rmdir \$LOCALDIR
"""
