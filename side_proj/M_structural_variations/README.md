## Structural variations

Detection of structural variants from population data will be based on pair-end illumina data.
One individual per species will have mate-pairs as well to access structural variations
with greater resolution.
Four more individuals of one species (_T. genevieve_) have PacBio reads, but we wont have more.

### Inference from Illumina reads

Methods are based on coverage or mapping patterns. Different methods produce different variant calls. This observation let to development of methods that calls consensus based on different methods. The consensus caller is called (Survivor)[https://github.com/fritzsedlazeck/SURVIVOR].

#### pair-end

##### Manta

The tools is designed for tumor cells, but there is no reason to think that it wont work on whole bodies.

I downloaded precompiled binaries of (Manta)[https://github.com/Illumina/manta] v1.0.3.

Unpacking and running `runMantaWorkflowDemo.py` reported succesful run and some details of commends that were tested. I will extract those commands and apply them to Tge, since I already have indexed reference and mapped pair end reads..

```
/Home/kjaron/src/manta-1.0.3.centos5_x86_64/bin/configManta.py --normalBam="map_pe_to_5_Tge.bam"  --referenceFasta="Tge_abyss87_besst_GC_core_100k_filtered.fa" --runDir=$(pwd)
```

the workflow script was created, so I run it. It can nicely limit the memory.

```
python -E runWorkflow.py -m local -j 32 -g 120 1> tge_350_manta.log 2> tge_350_manta.err
```

53k SVs, 7k of good quality. Very small Indels included.

Turned into scripts:

The name of the reference individual's insert size 350: `ref_is350`

```bash
manta.sh 1_Tdi b3v04 ref_is350
```

##### Delly

installed it is running

```
delly call -g Tge_abyss87_besst_GC_core_100k_filtered.fa \
  map_pe_to_5_Tge.bam -o tge_350_delly.bcf
```

produced in 2 hours a `.bcf` file that was not parsed so far.

##### Lumpy

requires `samtools` or `sambamba`.

##### Brakedancer

weird problem with libraries.

#### mate-pairs

I have mate-pairs for only reference genomes

### Inference from PacBio reads - does not matter

I am testing NGM-LR -> Sniffles pipeline; installation of both those tools was trivial! They both have man pages and Fritz was friendly, when he gave a talk about it, I guess I can consult it with him after some trials.

I should run it on
```
/scratch/beegfs/monthly/kjaron/5_Tge/asm_abyss_besst_gc/Tge_abyss87_besst_GC.fasta
```

but I run it on non-gapfilled version

```
/scratch/beegfs/monthly/kjaron/5_Tge/BESST_mapping/Tge_abyss87_besst.fa
```

```sh
bash $TROOT/N_variant_calling/2_map_longreads_lsf.sh Tge_GECD_map /scratch/beegfs/monthly/kjaron/5_Tge/BESST_mapping/Tge_abyss87_besst.fa /scratch/beegfs/monthly/kjaron/timema_PacBio_reads/5_Tge/filtered_subreads.GECD.7smrt.fastq.gz
/usr/bin/time -f '%M %E %P' sniffles -m Tge_GECD_map.bam -s 3 -t 4 -v GECD_calls.vcf
```

Crucial parameters seems to be `-s 3` which says how much coverage we require for accepting a structural variant, since the coverage of both datasets is small (<5x, I set it so far to 3).

Sniffles takes ages, preraps because I inputed too many scaffolds, I will remap reads to core scaffolds. I copied the result it computed so far, I will compare it to the fresh one.

```
cd /scratch/beegfs/monthly/kjaron/5_Tge/variant_calling
mkdir GEEF
cd GEEF
bash $TROOT/N_variant_calling/2_map_longreads_lsf.sh Tge_GEEF_map $TROOT/data/5_Tge/reference/Tge_abyss87_besst_GC_core_scaffolds.fa $TROOT/data/5_Tge/long_reads/filtered_subreads.GEEF.10smrt.fastq.gz
mkdir ../GECD
cd ../GECD
bash $TROOT/N_variant_calling/2_map_longreads_lsf.sh Tge_GECD_map $TROOT/data/5_Tge/reference/Tge_abyss87_besst_GC_core_scaffolds.fa $TROOT/data/5_Tge/long_reads/filtered_subreads.GECD.7smrt.fastq.gz
```

I still have a lot of scaffolds for Sniffles. I will try to make a subset of a few long scaffolds.

One job finished, but 55% of reads are not mapping to the reference. Second job have crushed in the middle, issue reported and job re-lunched with even more turnicated reference.

```
cd ..
mkdir Tge_GECD_remap
cd Tge_GECD_remap
python3 $TROOT/scripts/generic_genomics/fasta2fasta_length_filtering.py $TROOT/data/5_Tge/reference/Tge_abyss87_besst_GC_core_scaffolds.fa 50000 > Tge_abyss87_besst_GC_50k_filter.fa
bash $TROOT/N_variant_calling/2_map_longreads_lsf.sh Tge_GECD_remap $(pwd)/Tge_abyss87_besst_GC_50k_filter.fa $TROOT/data/5_Tge/long_reads/filtered_subreads.GEEF.10smrt.fastq.gz
```

and I started sniffles for GEEF

```
bash $TROOT/N_variant_calling/3_sniffles_lsf.sh Tge_GEEF_sv /scratch/beegfs/monthly/kjaron/5_Tge/variant_calling/GEEF/Tge_GEEF_map.bam
```

and GECD as well:

```
bash $TROOT/N_variant_calling/3_sniffles_lsf.sh Tge_GECD_sv $TROOT/data/5_Tge/variant_calling/remap_GECD/Tge_GECD_remap.bam
```
