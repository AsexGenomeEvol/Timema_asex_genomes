#/bin/bash

# 1 a directory with input reads
# 2 output dir data/{sp}/genomescope/{sample}

FASTQS=$(ls $1)
OUTDIR=$2
mkdir -p $OUTDIR

#zcat $1 $2 > $OUTDIR/trimmed_reads.fasta

jellyfish count -C -m 21 -s 1000000000 -t 16 -o $OUTDIR/kmer_counts <(zcat $FASTQS)

if [[ $(echo $OUTDIR/kmer_counts_* | wc -w) -eq 1 ]]
then
    mv $OUTDIR/kmer_counts_0 $OUTDIR/kmer_counts.jf
else
    jellyfish merge $OUTDIR/kmer_counts_* -o $OUTDIR/kmer_counts.jf && \
    rm $OUTDIR/kmer_counts_*
fi

jellyfish histo -t 16 $OUTDIR/kmer_counts.jf > $OUTDIR/kmer.hist
genomescope.R $OUTDIR/kmer.hist 21 100 $OUTDIR 1000 verbose

if [[ -s $OUTDIR/summary.txt ]]
then
    rm $OUTDIR/kmer_counts.jf
fi