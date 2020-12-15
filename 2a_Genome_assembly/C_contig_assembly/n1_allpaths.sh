#!/bin/bash


# default input files in_groups.csv in_libs.csv has to be present 
PrepareAllPathsInputs.pl DATA_DIR=`pwd`/Allpaths/1_Tps/DATA PLOIDY=2

# run allpaths-lg
RunAllPathsLG PRE=Allpaths REFERENCE_NAME=1_Tps DATA_SUBDIR=DATA RUN=default THREADS=32 1> allpaths.out 2> allpaths.err

# when it crushed...

RunAllPathsLG PRE=Allpaths REFERENCE_NAME=1_Tps DATA_SUBDIR=DATA RUN=default THREADS=32 OVERWRITE=True 1> allpaths_l2.out 2> allpaths_l2.err

# crushed for second time :-( seems that my overlapping library is not overlapping enough
