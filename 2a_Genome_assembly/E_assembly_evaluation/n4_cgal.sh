#!/bin/bash

#first argument: the assembly
ASMBLY=$1

#second argument: the reads mapped by bowtie2 to the assembly (in sam format)
MAP=$2

# IMPORTANT NOTE:
# cgal is saving files with names like "myout.sam" and every next step is expecting to find files with those names in the folder, where is it run.

# software is installed locally in /home/kjaron/bin
export PATH=$PATH:/home/kjaron/bin

# filter out aligments without * cigar strings tagged by YT:Z:UP (not-recognised by cgal for some unknown reason)
grep -v "YT:Z:UP" $MAP > filt_map.sam

# convert sam to cgal friendly sam
bowtie2convert filt_map.sam

# you can align unaligned reads, but it is very very slow (one day took to align 500 reads)
# align $ASMBLY <number of reads to align> <number of threadsa>

# finally compute likelihood
cgal $ASMBLY
