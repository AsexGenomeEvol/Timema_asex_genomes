#!/bin/bash

INBAM=$2
OUTBAM=$4

scripts/filter_splitreads.py $INBAM $OUTBAM
