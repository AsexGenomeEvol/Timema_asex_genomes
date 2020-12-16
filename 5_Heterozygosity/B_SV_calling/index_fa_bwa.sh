#!/bin/bash
#
# index fasta specified by 1st argument

module add UHTS/Aligner/bwa/0.7.17

bwa index $1
