#!/bin/bash
#
# is SP
# is a reference genome (VERSION)
# individual

module add UHTS/Analysis/delly/0.7.8

GENOME="$1"
BAM="$2"
DELLY="$4"

delly call -g $GENOME $BAM -o $DELLY