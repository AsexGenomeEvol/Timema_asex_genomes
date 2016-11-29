# Variant calling

I am interested in SNPs, small structural variants (based on illumina data) and greater structural variants (based on PacBio data).

## small

## big

I am testing NGM-LR -> Sniffles pipeline; installation of both those tools was trivial! They both have man pages and Fritz was friendly, when he gave a talk about it, therefore I guess I can consult it with him after several trials.

```sh
bash ~/timema_assembly/N_variant_calling/2_map_longreads_lsf.sh Tge_GECD_map /scratch/beegfs/monthly/kjaron/5_Tge/BESST_mapping/Tge_abyss87_besst.fa /scratch/beegfs/monthly/kjaron/timema_PacBio_reads/5_Tge/filtered_subreads.GECD.7smrt.fastq.gz
/usr/bin/time -f '%M %E %P' sniffles -m Tge_GECD_map.bam -s 3 -t 4 -v GECD_calls.vcf
```

Crucial parameters seems to be `-s 3` which says how much coverage we require for accepting a structural variant, since the coverage of both datasets is very small (<5x, I set it so far to 3).

