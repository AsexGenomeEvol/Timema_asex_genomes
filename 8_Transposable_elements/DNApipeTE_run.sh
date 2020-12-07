#!/bin/bash
#script to run dnaPipeTE

#get dnaPipeTE from https://github.com/clemgoub/dnaPipeTE
#follow installation instructions


#set up environment
#module add Development/java/latest
module add Development/java_jre/1.7.0_17
module add UHTS/Aligner/bowtie2/2.3.0
module add R/latest

export TMPDIR=/scratch/local/jbast/tmp/
export _JAVA_OPTIONS="-XX:ParallelGCThreads=24"

clear

#run the dnaPipeTE pipeline
#IMPORTANT: must be run from installed software path !!

cd /scratch/beegfs/monthly/jbast/software/dnaPipeTE_old/1.2/

#run with SE reads (can be gzipped)
#IMPORTANT: give genome size of organism

python3 ./dnaPipeTE.py -input /scratch/local/jbast/mites/On/On6.trim.sorted.filtered.180.pair1.fastq.gz \
 -output /scratch/local/jbast/dnaPipeTE/On \
-genome_size 193261443 -genome_coverage 0.5 -sample_number 4


module rm Development/java_jre/latest
module rm UHTS/Aligner/bowtie2/2.3.0
module rm R/latest

printf "****DONE****\n"
