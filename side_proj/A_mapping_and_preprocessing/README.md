

A code that would download raw reds of one resequencing ind from our ftp server. To get all of them, just wrap it around one more loop. Something like `for smpl in $SAMPLE; do ... done`.

```{bash}
for sp in $TIMEMAS; do
    SAMPLES=$(grep $sp data/resequencing_samples | cut -f 1 -d " ")
    smpl=$(echo $SAMPLES | cut -f 1 -d ' ')
    ODIR=data/"$sp"/raw_reads/$(grep "$smpl" data/resequencing_samples | cut -f 3)
    mkdir -p $ODIR
    wget ftp://ftpmrr.unil.ch/AsexGenomeEvol/timema/data/"$sp"/raw_reads/"$smpl"/* -P $ODIR
done
```

Trimming should have an extra level of complexity, something that takes species and sample name and trim all the reads.

