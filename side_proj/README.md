# Variant calling

I am interested in SNPs, small structural variants (based on illumina data) and greater structural variants (based on PacBio data).

#### SNPs & Indels

so far I was thinking of FreeBayes for variant calling, but according the latest benchmarking I should probably go for Octopus.

But quick solution is GATK for now. Mapping using `bwa-mem`

```
module add UHTS/Aligner/bwa/0.7.13

cd $TROOT/data/5_Tge/reference

# only scaffolds longer than 100k are kept (to make it computable)
python3 $TROOT/scripts/generic_genomics/fasta2fasta_length_filtering.py Tge_abyss87_besst_GC_core_scaffolds.fa 100000 > Tge_abyss87_besst_GC_core_100k_filtered.fa

# the total size of this reference is ~600M.
bwa index Tge_abyss87_besst_GC_core_100k_filtered.fa

$TROOT/scripts/map_pair_end_RG_lsf.sh 5_Tge \
  Tge_abyss87_besst_GC_core_100k_filtered.fa \
  Tge_R1t_is_350.fq.gz \
  Tge_R2t_is_350.fq.gz \
  "@RG\tID:TGE\tSM:TGE_REF\tPL:illumina\tLB:350\tPU:lane1"
```

and similar for `1_Tps`:

```
cd $TROOT/data/1_Tps/raw_assembly

# only scaffolds longer than 100k are kept (to make it computable)
python3 $TROOT/scripts/generic_genomics/fasta2fasta_length_filtering.py 1_Tps_genome.fa 100000 > 1_Tps_genome_100k_filtered.fa

# the total size of this reference is ~200M, 1326 scaffolds.
bwa index 1_Tps_genome_100k_filtered.fa

$TROOT/scripts/map_pair_end_RG_lsf.sh 1_Tps \
  Tge_abyss87_besst_GC_core_100k_filtered.fa \
  Tps_R1t_is_350.fq.gz \
  Tps_R2t_is_350.fq.gz \
  "@RG\tID:TPS\tSM:TPS_REF\tPL:illumina\tLB:350\tPU:lane2"
```

I am supposed to mark duplicates, but I wont, seems like a pint I do not have time for it now. However, script `1_deduplicate_lsf.sh` should do it (not tested).

Now we have to index fasta
```
cd $TROOT/data/1_Tps/reference
$TROOT/scripts/index_fa.sh 1_Tps_genome_100k_filtered.fa
$TROOT/scripts/make_dict_fasta.sh 1_Tps_genome_100k_filtered.fa
cd ../variant_calling/GATK/350/
$TROOT/scripts/index_bam.sh map_pe_to_1_Tps.bam
# once those jobs are done
$TROOT/N_variant_calling/2_run_GATK_lsf.sh 1_Tps 350 1_Tps_genome_100k_filtered.fa
```
analogically I run the analysis on Tge, both 350 and 550 libs.


#### Structural variations

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
