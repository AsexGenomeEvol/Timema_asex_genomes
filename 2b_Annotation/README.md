## Pipeline description for protein-coding genes annotation.

## Table of contents

1. [Required files](#1_required)
	* [Final genome](#11_genome)
	* [Busco](#12_busco)
	* [Transcript evidence](#13_transcriptome)	
	* [Protein evidence](#14_protein)

2. [Structural annotation](#2_structural)

3. [Functional annotation](#3_functional)

4. [EMBL format](#4_embl)


## <a name="1_required"></a>1) Required files 

#### <a name="11_genome"></a>1.1) Final genome

- `<sp>_b3v06` have been used.

#### <a name="12_busco"></a>1.2) BUSCO

[BUSCO](http://busco.ezlab.org/) has been run with the database `insecta_odb9`.
```
source Software/busco/3.0.2b.sh

run_BUSCO.py --long -i *.fasta -o output -l insecta_odb9 -m geno -c 10
```

#### <a name="13_transcriptome"></a>1.3) Transcript evidence

For the transcriptome: [STAR](https://github.com/alexdobin/STAR) for mapping and [Trinity genome guided](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Genome-Guided-Trinity-Transcriptome-Assembly) for the assembly. 

1) Map the reads with STAR (single end mode).

```
module add UHTS/Aligner/STAR/2.5.3a;

# Index creation
 
mkdir index
STAR --runMode genomeGenerate --genomeDir index --genomeFastaFiles genome.fasta --runThreadN 15

# Mapping step

mkdir 2pass_basic

STAR --genomeDir ../index/ --readFilesIn <*.fastq.gz ...> --runThreadN 10 --outSAMtype BAM SortedByCoordinate --twopassMode Basic --readFilesCommand zcat --limitBAMsortRAM 8051268437 --outFileNamePrefix 2pass_basic/out.
```

2) Assembly with Trinity (genome guided mode).

```
module add UHTS/Assembler/trinityrnaseq/2.4.0;

# out.Aligned.sortedByCoord.out.bam output by STAR.

Trinity --genome_guided_bam out.Aligned.sortedByCoord.out.bam --genome_guided_max_intron 100000 --max_memory 200G --CPU 35 --SS_lib_type R
```

Transcriptome unfiltered: `*Tv01.fasta`.

3) Filter the transcriptome.

a) Map all samples against transcriptome and compute the TPM value per sample using [Kallisto](https://pachterlab.github.io/kallisto/about). 

```
module add UHTS/Analysis/kallisto/0.43.0

kallisto quant -i <index_transcriptome> -t 10 --single -l 200 -s 20 --rf-stranded --bias -o <output> <reads>
```

b) Keep contigs that have at least 1 TPM in any sample.

```
module add R/3.3.2;

/software/UHTS/Assembler/trinityrnaseq/2.4.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix merge --name_sample_by_basedir <*abundance.tsv ...> 

/software/UHTS/Assembler/trinityrnaseq/2.4.0/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl merge.TPM.not_cross_norm | tee merge.TPM.not_cross_norm.counts_by_min_TPM

/software/UHTS/Assembler/trinityrnaseq/2.4.0/util/filter_low_expr_transcripts.pl --matrix merge.TPM.not_cross_norm --transcripts *Tv01.fasta --min_expr_any 1 > *Tv02.fasta
```

Transcriptome filtered: `*Tv02.fasta`.

#### <a name="14_protein"></a>1.4) Protein evidence

- UniProtKB/Swiss-Prot (release 2018_01) + BUSCO `insecta_odb9`.

```
cat uniprot_sprot.fasta insecta_odb9/ancestral_variants > uniprot_sprot_insecta.fasta
```

## <a name="2_structural"></a>2) Structural annotation 

[MAKER](http://www.yandell-lab.org/software/maker.html) (v. 2.31.8) has been used.

Set the environment (Only for Vital-IT):

```
export PATH=/scratch/beegfs/monthly/ptranvan/Software/mpich-3.2/mpi_install/bin:/scratch/beegfs/monthly/ptranvan/Software/maker_7_local/bin:/software/SequenceAnalysis/GenePrediction/snoscan/0.9.1/bin:/software/SequenceAnalysis/GenePrediction/tRNAscan-SE/2.0.0/bin:/software/SequenceAnalysis/GenePrediction/augustus/3.2.3/scripts:/software/SequenceAnalysis/GenePrediction/augustus/3.2.3/bin:/scratch/beegfs/monthly/ptranvan/Software/RepeatMasker/4.0.7_local:/software/Blast/ncbi-blast/2.7.1+/bin:/software/SequenceAnalysis/Repeat/trf/4.07b/bin:/software/SequenceAnalysis/HMM-Profile/hmmer/3.1b2/bin:/software/SequenceAnalysis/Repeat/RMBlast/2.6.0+/bin:/software/UHTS/Aligner/phrap/0.990329/bin:/software/SequenceAnalysis/SequenceAlignment/exonerate/2.4.0/bin:/software/SequenceAnalysis/GenePrediction/snap/2013.11.29/bin:/home/ptranvan/qiime/bin/:/home/ptranvan/qiime/bin/:/mnt/common/lsf/9.1/linux2.6-glibc2.3-x86_64/etc:/mnt/common/lsf/9.1/linux2.6-glibc2.3-x86_64/bin:/software/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/scratch/beegfs/monthly/ptranvan/Software/assemblathon/1.0:/home/ptranvan/bin:/scratch/beegfs/monthly/ptranvan/Software/assemblathon/1.0

export AUGUSTUS_CONFIG_PATH=/scratch/beegfs/monthly/ptranvan/Software/busco/3.0.2b/augustus_config_7
```

#### <a name="21_round1"></a>2.1) Iteration 1


1) Create and modify control files. See [config_files](./config_files/) for more details.

```
mkdir maker
cd maker

maker -CTL
```

Parameters to change:

- In `maker_bopts.ctl`:

```
depth_blastn=10
depth_blastx=10
depth_tblastx=10
```

- In `maker_opts.ctl`:

```
genome=../*b3v06.fasta
est=../evidence/*_Tv02.fasta
protein=../evidence/uniprot_sprot_insecta.fasta

augustus_species=BUSCO_*
est2genome=1
protein2genome=1

max_dna_len=300000
alt_splice=1
split_hit=100000
correct_est_fusion=1
```
 
3) Run MAKER.

```
mpiexec -n 20 maker
```

Once finished:

```
gff3_merge -d *.maker.output/*_master_datastore_index.log

mv maker_opts.ctl maker_opts_round1.ctl
mv *.all.gff *.all.round1.gff
cp maker_opts_round1.ctl maker_opts.ctl
```

#### <a name="22_round2"></a>2.2) Iteration 2


1) Train MAKER result from iteration 1 with [SNAP](https://github.com/KorfLab/SNAP). 

```
mkdir -p snap/round1
cd snap/round1

maker2zff ../../maker/*.all.round1.gff 

fathom genome.ann genome.dna -categorize 1000
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna

hmm-assembler.pl * . > *.maker_round1_snap.hmm
```

2) Train MAKER result from iteration 1 with [Augustus](http://augustus.gobics.de/).

```
autoAug.pl --genome=../snap/round1/export.dna --species=AUGUSTUS_*_round1 --trainingset=../snap/round1/export.aa --singleCPU -v --useexisting > output.log
```

3) Modify control files.

Parameters to change:

- In `maker_opts.ctl`:

```
snaphmm=../snap/round1/*.maker_round1_snap.hmm
augustus_species=AUGUSTUS_*_round1
est2genome=0
protein2genome=0
```
 
3) Run MAKER.

```
mpiexec -n 20 maker
```

Once finished:

```
gff3_merge -d *.maker.output/*_master_datastore_index.log

mv maker_opts.ctl maker_opts_round2.ctl
mv *.all.gff *.all.round2.gff
cp maker_opts_round2.ctl maker_opts.ctl
```

## <a name="3_functional"></a>3) Functional annotation 


1) Create GFF file.
 
```  
mkdir functional
cd functional

awk '/\tmaker\t/' ../*.all.round2.gff > *.max.gff
   
```    
    
2) Create protein and transcript fasta files.
  
```    
fasta_merge -d ../*.maker.output/*_master_datastore_index.log

cp *.all.maker.proteins.fasta *.max.proteins.fasta
cp *.all.maker.transcripts.fasta *.max.transcripts.fasta
``` 

3) Assign short id to each protein coding genes. 

``` 
maker_map_ids --prefix TCM_ --justify 5 *.max.gff > *.max.map	#example for T. californicum.
    
map_gff_ids *.max.map *.max.gff
map_fasta_ids *.max.map *.max.proteins.fasta
map_fasta_ids *.max.map *.max.transcripts.fasta
``` 

4) Assign gene functions and GO-terms with [Blast2GO](https://www.blast2go.com/).

Default parameters against both the `NCBI non-redundant arthropods protein` database and the Drosophila melanogaster `drosoph` database. 

## <a name="4_embl"></a>4) EMBL format 

The genome and the annotation have been converted to EMBL format in order to be deposited in the European Nucleotide Archive (ENA) using [AGAT](https://github.com/NBISweden/AGAT) and [EMBLmyGFF3](https://github.com/NBISweden/EMBLmyGFF3).

1) Fixing mRNA duplicates and flag genes with short intron.

```
agat_sp_fix_features_locations_duplicated.pl --gff * --out *

agat_sp_flag_short_introns.pl --gff * --out *
```

2) Convert to EMBL format.

```
EMBLmyGFF3 *.gff *.fasta --topology linear --molecule_type 'genomic DNA' --transl_table 1  --species '*' --locus_tag * --project_id * --de * -o *
```



