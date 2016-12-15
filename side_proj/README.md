# Variant calling

I am interested in SNPs, small structural variants (based on illumina data) and greater structural variants (based on PacBio data).

## small

## Variant calling

so far I was thinking of FreeBayes for vartiant calling, but according the latest benchmarking I should probably go for Octopus (which is not available now).


## big

I am testing NGM-LR -> Sniffles pipeline; installation of both those tools was trivial! They both have man pages and Fritz was friendly, when he gave a talk about it, therefore I guess I can consult it with him after several trials.

I should run it on
```
/scratch/beegfs/monthly/kjaron/5_Tge/asm_abyss_besst_gc/Tge_abyss87_besst_GC.fasta
```

but unfortunately I accidentally run it on non-gapfilled version

```
/scratch/beegfs/monthly/kjaron/5_Tge/BESST_mapping/Tge_abyss87_besst.fa
```

```sh
bash ~/timema_assembly/N_variant_calling/2_map_longreads_lsf.sh Tge_GECD_map /scratch/beegfs/monthly/kjaron/5_Tge/BESST_mapping/Tge_abyss87_besst.fa /scratch/beegfs/monthly/kjaron/timema_PacBio_reads/5_Tge/filtered_subreads.GECD.7smrt.fastq.gz
/usr/bin/time -f '%M %E %P' sniffles -m Tge_GECD_map.bam -s 3 -t 4 -v GECD_calls.vcf
```

Crucial parameters seems to be `-s 3` which says how much coverage we require for accepting a structural variant, since the coverage of both datasets is very small (<5x, I set it so far to 3).

Sniffles takes ages, progrably because I inputed too many scaffolds, I will remap reads to just core scaffolds. I copied the result it computed so far, I will compare it to the fresh one.

```
cd /scratch/beegfs/monthly/kjaron/5_Tge/variant_calling
mkdir GEEF
cd GEEF
bash ~/timema_assembly/N_variant_calling/2_map_longreads_lsf.sh Tge_GEEF_map /scratch/beegfs/monthly/kjaron/5_Tge/asm_abyss_besst_gc/Tge_abyss87_besst_GC_core_scaffolds.fa /scratch/beegfs/monthly/kjaron/timema_PacBio_reads/5_Tge/filtered_subreads.GEEF.10smrt.fastq.gz
mkdir ../GECD
cd ../GECD
bash ~/timema_assembly/N_variant_calling/2_map_longreads_lsf.sh Tge_GECD_map /scratch/beegfs/monthly/kjaron/5_Tge/asm_abyss_besst_gc/Tge_abyss87_besst_GC_core_scaffolds.fa /scratch/beegfs/monthly/kjaron/timema_PacBio_reads/5_Tge/filtered_subreads.GECD.7smrt.fastq.gz
```

but it seems that I still have a lot lot scaffolds for Sniffles. I will try to make a subset of really few very long scaffolds.

One job finished, seems successfully, but 55% of reads are not mapping to the reference. Second job have crushed in the middle, issue reported and job relunched with even more turnicated reference.

```
cd ..
mkdir Tge_GECD_remap
cd Tge_GECD_remap
python3 ~/scripts/generic_genomics/fasta2fasta_length_filtering.py /scratch/beegfs/monthly/kjaron/5_Tge/asm_abyss_besst_gc/Tge_abyss87_besst_GC_core_scaffolds.fa 50000 > Tge_abyss87_besst_GC_50k_filter.fa
bash ~/timema_assembly/N_variant_calling/2_map_longreads_lsf.sh Tge_GECD_remap $(pwd)/Tge_abyss87_besst_GC_50k_filter.fa /scratch/beegfs/monthly/kjaron/timema_PacBio_reads/5_Tge/filtered_subreads.GEEF.10smrt.fastq.gz
```

and I started sniffles for GEEF

```
bash ~/timema_assembly/N_variant_calling/3_sniffles_lsf.sh Tge_GEEF_sv /scratch/beegfs/monthly/kjaron/5_Tge/variant_calling/GEEF/Tge_GEEF_map.bam 
```
