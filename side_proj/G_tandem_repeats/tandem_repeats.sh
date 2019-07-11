#!/bin/bash

# all but last argument
READS=${@:1:$#-1}
# last argument
OUTFILE=${@:$#}
# remove suffix of the last argument
OUTPATTERN=${OUTFILE%.*}

zcat $READS > temp.fastq

perl /home/kjaron/bin/k_seek.pl temp.fastq $OUTPATTERN

rm temp.fastq

echo "Done"