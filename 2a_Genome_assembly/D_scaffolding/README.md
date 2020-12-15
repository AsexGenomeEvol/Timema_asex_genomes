## Testing and details

Here I show test using only mate pairs, but in fact in the latest version I have used also pair-end reads. Using also pair-end reads improved scaffolding over tests.

#### BESST on ABySS

is a software that handles pair-end contamination in mate-pair library. That allow
to be more relaxed about mate-pairs filtering.

Tests were done at commit [981a20](https://github.com/AsexGenomeEvol/timema_assembly/tree/981a203e91ef5ecab82950369cee452b53e8d5b2/D_scaffolding).
BESST script for mapping mate-pairs using bwa-mem on `5_Tge` contigs.

```bash
module add UHTS/Aligner/bwa/0.7.13
/usr/bin/time -f '%M %E %P' python ~/src/BESST/scripts/reads_to_ctg_map.py Tge_R1t_is_3000.fq.gz Tge_R2t_is_3000.fq.gz ``/scratch/beegfs/monthly/kjaron/5_Tge/abyss_k83/Tge_abyss_k83-6.fa --threads 16
```

this assembly was scaffolded and gapfilled using GapCloser and than filtered using [the script](../C_contig_assembly/coverage_filtering.R). The core 900M in 168k scaffolds had almost all busco genes using [this script](../E_assembly_evaluation/2_busco.sh):

```
Complete:91.7% [Single: 91.3%, Duplicated:0.4%]
Fragmented:5.3%,
Missing:3.0%,
number of busco genes:1658
```

on `2_Tcm`

```sh
$TROOT/scripts/map_mate_pairs_lsf.sh 2_Tcm_map_is3000 `pwd`/Tcm_abyss_k83_kc2_filt200.fa \
  $TROOT/data/2_Tcm/trimmed_reads/mp_nxtrim_RF/Tcm_R1t_is_3000.fq.gz \
  $TROOT/data/2_Tcm/trimmed_reads/mp_nxtrim_RF/Tcm_R2t_is_3000.fq.gz
$TROOT/scripts/map_mate_pairs_lsf.sh 2_Tcm_map_is5000 `pwd`/Tcm_abyss_k83_kc2_filt200.fa \
  $TROOT/data/2_Tcm/trimmed_reads/mp_nxtrim_RF/Tcm_R1t_is_5000.fq.gz \
  $TROOT/data/2_Tcm/trimmed_reads/mp_nxtrim_RF/Tcm_R2t_is_5000.fq.gz
$TROOT/D_scaffolding/1_BESST_lsf.sh 2_Tcm_k83_BESST \
  `pwd`/Tcm_abyss_k83_kc2_filt200.fa `pwd`/2_Tcm_map_is3000.bam `pwd`/2_Tcm_map_is5000.bam
```

the script has some issues. To debug it:

```bash
cd $TROOT/data/2_Tcm/assembly/abyss_k83_correctedreads
python ~/src/BESST/runBESST -c 2_Tcm_abyss_k83-6.fa -f 2_Tcm_map_is3000.bam 2_Tcm_map_is5000.bam -o 2_Tcm_k83cr_BESST.fa --orientation rf rf -m 3000 5500 -s 450 1400
```

on the same data to get better scaffolds:

```bash
cd /scratch/beegfs/monthly/kjaron/timema_assembly/data/2_Tcm/assembly/scaffolding
python ~/src/BESST/runBESST -c Tcm_abyss_k83_kc2_filt200.fa -f 2_Tcm_map_is3000.bam 2_Tcm_map_is5000.bam -o 2_Tcm_k83_BESST_run2 -orientation rf rf -K 83 -plots -filter_contigs 500
```

ok, this one was total disaster. I guess I should filter the contig assembly first. I know already that it is substantially longer than it should be and the most of those contigs are of a size of a read. I will just assume that the most of the information got in contigs longer and 200 bases, that will make life of the scaffolder way way easier! I can think of a way how to remove duplicates from assembly afterwards (in case that after scaffolding the length will be too high)

`ABySS` can be also used for estimation of real insert sizes of mate pairs (can be found in `.err` log).

#### SOAP scaffolding of ABySS assembly

ok there is this really strange thing going on - ABySS has reasonable contigs and relatively fragmented scaffolds (N50 ~ 18k) and SOAP has fragmented contigs (only hald of the genome is covered by contigs > 500) and substantially more continuity (N50 ~ 35k). What happens when I scaffold abyss contigs with SOAP scaffolder??

```bash
SOAPdenovo-63mer map -s 2_Tcm.config -g Tcm_abyss_SOAP_graph 1>SOAP_map.log 2>SOAP_map.err
```

nope, it is not working. At least this version requires index of contigs that seems not to be stadardised and at the same time they do not provide a way now to create it. However, we have this several years old version on cluster.

Ha, v251 offers a new binary: `SOAPdenovo-fusion` that suppose to allow it. I will try to scaffold it with high kmer first.

```bash
cd /scratch/local/kjaron/SOAP_scfder
./SOAPdenovo-fusion -D -K 83 -c Tcm_abyss_k83_kc2_filt200.fa -g 2_Tcm_scf_k83 1> SOAPfusion.out 2> SOAPfusion.err
./SOAPdenovo-127mer map -s 2_Tcm.config -g 2_Tcm_scf_k83 -p 32 1> SOAP_map.log 2> SOAP_map.err
./SOAPdenovo-127mer scaff -g 2_Tcm_scf_k83 -p 32 1> SOAP_scaff.log 2> SOAP_scaff.err
```

the result was very very bad (fragmented).

##### SOAP scaffolding of Megahit assembly

```bash
cd /scratch/local/kjaron/SOAP_scfder
export TMPDIR=`pwd`/temp
ln -s /scratch/local/kjaron/2_Tcm_corrected_reads/Tcm_megahit_k47-83_run1/megahit_out/Tcm_megahit_k47-83.contigs.fa .
./SOAPdenovo-fusion -D -K 47 -c Tcm_megahit_k47-83.contigs.fa -g 2_Tcm_megahit_SOAP_scf_k47 &> 2_Tcm_megahit_SOAP_fusion_k47.err
./SOAPdenovo-127mer map -s 2_Tcm.config -g 2_Tcm_megahit_SOAP_scf_k47 -p 32 &> 2_Tcm_megahit_SOAP_map_k47.err
./SOAPdenovo-127mer scaff -g 2_Tcm_megahit_SOAP_scf_k47 -p 32 -N 1300000000 &> 2_Tcm_megahit_SOAP_scaf_k47.err
```

and very same for k = 83. The multi kmer approach yielded with assemblies inferior to single kmer SOAPdenovo2 pipeline.

#### Scaffold filtering

lots of contigs have very strange coverages - I have three libraries, which I use to filter out all very small contigs and contigs with high variance in coverage.
Remaining contigs are categorised to two groups - **funky** (where coverage in at least two libraries is lower than 0.5 or higher than 1.5 of the median computed from the coverage of contigs > 1k). The rest of the contigs form a second group called **core** cotigs.

Tests are run on Tge - scaffolds (yes, I think it would be better to do it before scaffolding).

```bash
bash $TROOT/scripts/map_pair_end_lsf.sh Tge_pe_350_map $TROOT/PATHTO/Tge_abyss87_besst_GC.fasta $TROOT/data/5_Tge/timema_trimmed/Tge_R1t_is_350.fq.gz $TROOT/data/5_Tge/timema_trimmed/Tge_R2t_is_350.fq.gz
bash $TROOT/scripts/map_pair_end_lsf.sh Tge_pe_550 $TROOT/PATHTO/Tge_abyss87_besst_GC.fasta $TROOT/data/5_Tge/timema_trimmed/Tge_R1t_is_550.fq.gz $TROOT/data/5_Tge/timema_trimmed/Tge_R2t_is_550.fq.gz
bash $TROOT/scripts/map_pair_end_lsf.sh Tge_pe_700 $TROOT/PATHTO/Tge_abyss87_besst_GC.fasta $TROOT/data/5_Tge/timema_trimmed/Tge_R1t_is_700.fq.gz $TROOT/data/5_Tge/timema_trimmed/Tge_R2t_is_700.fq.gz
```

then locally I computed per contig coverage using `samtools` and [awk script](https://github.com/KamilSJaron/generic_genomics/blob/master/depth2depth_per_contig.awk)

```bash
samtools depth Tge_pe_350_map.bam | ~/scripts/generic_genomics/depth2depth_per_contig.awk > Tge_350_per_contig_depth.tsv
samtools depth Tge_pe_550.bam | ~/scripts/generic_genomics/depth2depth_per_contig.awk > Tge_550_per_contig_depth.tsv
samtools depth Tge_pe_700.bam | ~/scripts/generic_genomics/depth2depth_per_contig.awk > Tge_700_per_contig_depth.tsv
```

on the three output files I run an R script to create list of **funky** and **core** scaffolds.
Besides the list several plots documenting the filtering step are produced.

```bash
Rscript coverage_filtering.R
```

R script has to be bit generalized, so far it suited for Tge (and maybe it should be filtered in python instead for the performance reasons).

In the end I have not used this script to sort scaffolds, since it requires a lot of computational effort and it is very similar to filtering by a simple cut off 1000 bases.
