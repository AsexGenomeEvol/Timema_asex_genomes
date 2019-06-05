#!/bin/bash

GENOME="$1"
BAM="$4"
OUTDIR="$6"

# generate manta script
/Home/kjaron/src/manta-1.5.0.centos6_x86_64/bin/configManta.py \
  --normalBam="$BAM" \
  --referenceFasta="$GENOME" \
  --runDir="$OUTDIR"

python -E "$OUTDIR"/runWorkflow.py -m local -j 32 -g 180

# clean the temp data (tons and tons of files)
rm -r "$OUTDIR"/workspace/pyflow.data/logs/tmp
# tar the rest of the workspace
tar -cf "$OUTDIR"/workspace.tar "$OUTDIR"/workspace --remove-files

echo "Done"