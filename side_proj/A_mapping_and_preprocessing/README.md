

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

I am also adding samples `Tsp_00` that are libraries from the reference genome:

```
for sp in $TIMEMAS; do ln -s /scratch/beegfs/monthly/kjaron/timema_assembly/data/$sp/trimmed_reads/is_550 /scratch/beegfs/monthly/kjaron/variant_analysis/data/$sp/trimmed_reads/$(echo $sp | cut -f 2 -d _)_00; done
```