#!/bin/bash

# 1st argument if the species

if [[ ! -s $TROOT/data/$1/trimmed_reads/$1""_R1t_is_225.fq.gz ]]
then
	>&2 echo FILE $TROOT/data/$1/trimmed_reads/$1""_R1t_is_225.fq.gz IS MISSING;
	exit 1
fi

if [[ ! -s $TROOT/data/$1/trimmed_reads/$1""_R1t_is_350_run2.fq.gz ]]
then
        >&2 echo FILE $TROOT/data/$1/trimmed_reads/$1""_R1t_is_350_run2.fq.gz IS MISSING;
        exit 1
fi

if [[ ! -s $TROOT/data/$1/trimmed_reads/mp_nxtrim_FR/$1""_R1t_is_3000.fq.gz ]]
then
	>&2 echo FILE $TROOT/data/$1/trimmed_reads/mp_nxtrim_FR/$1""_R1t_is_3000.fq.gz IS MISSING;
	exit 1
fi

if [[ ! -s $TROOT/data/$1/trimmed_reads/mp_fasteris_FR/$1""_R1t_is_3000.fq.gz ]]
then
        >&2 echo FILE $TROOT/data/$1/trimmed_reads/mp_fasteris_FR/$1""_R1t_is_3000.fq.gz IS MISSING;
        exit 1
fi

if [[ $(ls /scratch/beegfs/monthly/kjaron/timema_assembly/data/$1""/trimmed_reads/*gz | wc -l) -ne 19 ]]
then
	>&2 echo ERROR: unxeptected number of trimmed pair-end files: $(ls /scratch/beegfs/monthly/kjaron/timema_assembly/data/$1""/trimmed_reads/*gz | wc -l) instead of 19;
	exit 1
fi

exit 0;
