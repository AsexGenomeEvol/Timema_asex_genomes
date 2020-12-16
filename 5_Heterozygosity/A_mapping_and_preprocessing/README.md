

A code that would download raw reds of one resequencing ind from our ftp server. To get all of them, just wrap it around one more loop. Something like `for smpl in $SAMPLE; do ... done`.

```{bash}
for sp in $TIMEMAS; do
    SAMPLES=$(grep $sp A_mapping_and_preprocessing/resequencing_samples | cut -f 1 -d " ")
    smpl=$(echo $SAMPLES | cut -f 1 -d ' ')
    ODIR=data/"$sp"/raw_reads/$(grep "$smpl" A_mapping_and_preprocessing/resequencing_samples | cut -f 3)
    mkdir -p $ODIR
    wget ftp://ftpmrr.unil.ch/AsexGenomeEvol/timema/data/"$sp"/raw_reads/"$smpl"/* -P $ODIR
done
```

Trimming was done in Montpellier, to sort out trimmed reads

```{bash}
for sp in $TIMEMAS; do
    sp_name=$(echo $sp | cut -f 2 -d "_");
    for order in 0{1,2,3,4,5}; do
        sample="$sp_name"_"$order";
        read_dir=data/"$sp"/trimmed_reads/"$sample";
        mkdir -p "$read_dir";
        mv data/reads_reseq_trimmed/"$sp"_"$order"_* $read_dir;
    done;
done

```

Alright, all the SV calls take for ever detecting breakpoint between scaffolds as SVs (scaffolds are considered chromosomes in the SV world). To prevent this unintended behaviour I need to get rid of all reads mapping to edges of different scaffolds.
I was considering `BAMQL`, but it was too complicated to install on the cluster.
So, I will map `A_mapping_and_preprocessing/filter_splitreads.py`
