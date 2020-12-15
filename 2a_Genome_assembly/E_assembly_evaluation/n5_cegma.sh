#!/bin/bash

#BSUB -L /bin/bash
#BSUB -q long
#BSUB -J cegma_baby
#BSUB -R "rusage[mem=5000]"
#BSUB -M 5000000


## running on vit_it


## cp genome to scratch - (change directory to your scratch!)
cp /archive/dee/schwander/dparker/timema_asm/*.fa /scratch/beegfs/monthly/dparker/

module load SequenceAnalysis/CEGMA/2.5

### runs cegma on all assemblies. Takes ages. Would be better to run all in parallel by adding a & to the end of the cegma -g "$foo1" -o "$foo2""_filt_300_cegma_out" line
for i in /scratch/beegfs/monthly/dparker/*.fa; do
    foo1=`echo "$i"`
    foo2=`echo "$i" | sed 's/.fa//'`
    cegma -g "$foo1" -o "$foo2""_filt_300_cegma_out"
done

