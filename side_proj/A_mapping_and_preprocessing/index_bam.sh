#!/bin/bash
#
# index bam specified by 1st argument

module add UHTS/Analysis/samtools/1.3

samtools index $1
