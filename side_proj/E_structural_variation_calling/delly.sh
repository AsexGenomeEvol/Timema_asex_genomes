#!/bin/bash
#
# is SP
# is a reference genome (VERSION)
# individual

GENOME="$1"
BAM="$4"
DELLY="$6"

# -r 200 is basically excluding translocations, that are impossible to detect due to heavy fragmentation of our references
delly call -r 200 -g $GENOME $BAM -o $DELLY

# we could possibly run the SV calls on Nosil's Tce reference, so we know how many translocations we are missing in other species or something...
