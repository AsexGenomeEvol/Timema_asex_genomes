# Assembly of Timema species

This is documentation of a procedure of genome assembly of the ten stick insect species. The references build in sections A-F has code `<sp>_b3v08`. All the downstream analysis are (going to be) based on it.

## Table of content:

TODO

## Organisation

The whole pipe line is in this README file, individual step has an individual subdirectory with scripts, and operates using following four directories:

- [data](data) - the big data (not present in this repository)
- [figures](figures) - figures
- [stats](stats) - small data, stats and overviews
- [scripts](scripts) - misc scripts

more technical documentation of each step is in every subfolder in their corresponding `README.md` files.

### Versioning of assemblies

We tested several assembly approaches (internally called "batches"), the only relevant one is based on ABySS and BESST pipeline (details bellow). The assembly went though several stages, every stage got assigned version name by the order.
The names of working versions of assemblies will contain information about

- `X` : sp pair [1 - 5]
- `sp` : species
- `ZZ`. version [00 - 99]

The format of name is `X_Tsp_b3vZZ.fa`. For example `1_Tdi_b3_v03` would be assembly of Timema douglasi (spaecies pair 1) from ABySS/BESST pipeline, version 03 (scaffolds).
Versioning of genomes will be described in following table with details in appropriate sections:

|  version  |    b3                     |
|:---------:| ------------------------- |
|    v01    |  contigs                  |
|    v02    |  filtered contigs         |
|    v03    |  scaffolds                |
|    v04    |  gap filled v03           |
|    v05    |  filtered v04 (> 1kbp)    |
|    v06    |  blob filtered v05        |
|    v07    |  blob filt v04 softmasked |
|    v08    |  v07 formatted for db     |

What the versions are for :
 - `v04` is for mapping reference reads, analysis is however performed only on a subset
 - `v06` is the subset of `v04` that is used for all the downstream analysis; it can serve as an input for analyses that do not involve read mapping
 - `v07` is the complete _Timema_ genome, but tiny scaffolds (<= 1000bp) are still potential contaminants.
 - `v08` is the complete _Timema_ genome, should be used for everything; this version is going to submitted to databases

### Environmental variable used in the workflow

To avoid absolute paths in scripts, I have decided to use environmental variable
TROOT (timema root) at the root folder of timema assembly project. Scripts are still
designed for cluster running `lsf` scheduler.

```bash
echo "export TROOT=/scratch/beegfs/monthly/kjaron/Timema_asex_genomes/2a_Genome_assembly" >> ~/.bashrc
```

General dependencies (used versions):
- `GNUmake` (3.81)
- `bash` (4.1.2(1)-release)
- `R` (3.3.2)
- `python` (3.4.1) with `Biopython`
- `bwa-mem` (0.7.15)
- `samtools` (1.3)

Dependencies by sections:
- [B_read_parsing](B_read_parsing) : `trimmomatic` (0.33), `NxTrim` (v0.4.1)
- [C_contig_assembly](C_contig_assembly) - `SOAPdenovo2` (2.04), `ABySS` (2.0.0)
- [D_scaffolding](D_scaffolding) - `SOAPdenovo2`, `ABySS`, `BESST` (2.2.6)
- [E_assembly_evaluation](E_assembly_evaluation) - `quast` (v4.1), `BUSCO`, `cutadapt` (v1.13), `Blobtools`

## Pipeline

**Notes about the organization of this pipeline**

- computed summary tables in stats directory are present in repository as well
- the make pipeline is semiautomatic. It means that you can't redo the whole study in one `make`, you need to build piece by piece using individual make instructions and you have to make sure that the next instruction is submitted once the previous one has successfully finished. This comes from two problems, first I have not chosen the optimal build language, GNU make is not aware of cluster computing, therefore it is non-trivial to write proper recipes (the choice was made at the time Snakemake was not so fun to work with). The second aspect is the size of the project, to redo everything top to bottom including all the tests you would need several computational years and sufficient resources. For this reason the default build target for make is its help, where all pipeline recipes are listed.
- Makefiles themselves are hard to read, if you are not sure what any step is doing, you can run `make <step> --dry-run`, which will just list command make would execute (shorter version of the same is `make <step> -n`), which is way easier than trying to fish from the makefile what it is supposed to do. Then you can look at the scripts written in friend languages like `R`, `python` or `bash`.

### Read parsing

The raw reads present `data/<sp>/raw_reads/<is>/<files>` need to be cleaned before the assembly. We aim to cut all low quality bases and and sequences of adapters using `trimmomatic`. Mate-pair libraries also require identification of linker sequence and reversing orientation to be Forward-Reverse direction, detailed information can be found at [webpage of Illumina](https://www.illumina.com/documents/products/datasheets/datasheet_genomic_sequence.pdf). We used `NxTrim` for delinking.

The trimmed pair-end reads are saved to

```
data/<sp>/trimmed_reads/<files>
```

and delinked and filtered mate-pair reads are saved to

```
data/<sp>/trimmed_reads/mp_nxtrim_FR/<files>
```

Every pair-end library is represented by four files, `*R1t / *R2t` - trimmed forward and reverse reads and `*R1np / *R2np` - single end reads surviving trimming without their pair read (single end cut from `R1` or `R2` respectively).
The additional sequencing of `is_350` library done eight months after the first sequencing round are saved with suffix `*_run2`.

The `is225` library is a library of pair-end contamination of mate-pair sequencing,
`se_mp` is a file gathering all single end reads surviving trimming without their mate-pair.

The optimal kmer for assembly of contigs was estimated using `kmergenie`.

#### execution

trim everything by (submits a lot of jobs to a cluster using `lsf`). Make is internally using scripts [trim_pair_end_reads_lsf.sh](B_read_parsing/trim_pair_end_reads_lsf.sh)
and [process_mate_pair_reads.sh](B_read_parsing/process_mate_pair_reads.sh)

```
make trimmed_reads
```

trimming will take a while. Once raw read processing is done you can clean raw reads, since they are not needed from now on:

```
make clean.raw_reads
```

To get an idea of the coverage after trimming we calculate number of nucleotides in every file

```
make stats.trimmed_reads
```

once all individual stats are computed, we pull counts into one [coverage table](stats/reads/coverage.table.tsv)

```
make stats/reads/coverage.table.tsv
```

Estimate optimal kmer for assembly

```
make kmergenie
```

The assembly of trimmed reads using optimal predicted kmer is described in directory [C_contig_assembly](../C_contig_assembly)

### Contig assembly

I wanted to separate steps of contig building and scaffolding to be able to optimize these two steps separately. However, I found that config assembly and scaffolding is more interlinked that I expected. For example `SOAPdenovo2` produces very fragmented contig assembly with size way over the size of genome, however SOAP scaffolder is apparently able to deal with it and in the end the scaffolds produced by `SOAPdenovo2` had continuity comparable to other pipelines.

I tested with more or less effort `Platanus`, `SOAPdenovo2`, `ABySS`, `Megahit`. I also tried `MaSuRCA` and `ALLPATHS-LG`, but I encountered technical difficulties. Beside a choice of assembler I explored effects of input data, namely the effect of correction of reads, effect of additional sequencing, adding single end reads, and using pair-end contamination of mate pairs. Based on trials I found three approaches that were tested on all ten species (see [C_contig_assembly](C_contig_assembly)). In the pipeline section only the final approach using `ABySS` is presented. The optimal kmer for `ABySS` was determined by `kmergenie` (and experimentally verified that it represent locally optimal assembly continuity).

#### execution

To assemble the genomes with `ABySS` run

```
make assemble.batch3
```

make stats of all the assemblies

```
make asm.stats
make asm.tables
```

Scaffolding of ABySS contigs and more tested scaffolding strategies are described in directory [D_scaffolding](D_scaffolding).

### Scaffolding and Gap closing

partially this section is covered in the previous section, I considered [BESST](https://github.com/ksahlin/BESST) and inbuilt scaffolders of all the individual assembly pipelines I tested (check README in the [D_scaffolding](D_scaffolding) directory for other tests). `BESST` has shown to over-perform all other approaches. The main advantage is that `BESST` considers that a certain proportion of mate pair libraries is just a pair-end contamination.

We also filled gaps of the scaffolds using `GapCloser` a part of SOAPdenovo package.

#### execution

`BESST` pipeline does not handle well too many contigs on input, therefore prior scaffolding we filter out everything smaller than 250 bases (couple of minutes). Then we build a bwa index (couple of hours) and finally map all pair end and mate pairs to contigs. Note that this is 5 (libraries) * 10 (species) = 50 jobs, where each job takes couple of hours. In the end we need to build an index for every single one of the bam files. Note that all this depends on existing files, so even the make files will not like that the files and directories do not exist.

```
make filter.batch3
make indexed.batch3
make mapping.batch3
make index.all.bam
```

scaffold `ABySS` contigs with `BESST` using mapped reads

```
make scaffold.batch3
```

At this level I see that batch2 (the one with more data) is actually bit worse than batch1, therefore I will not even try to gapfill it, gapifilling will be done only for batches 1 and 3.

```
make gapfilling
```

Detailed evaluation of scaffolded assemblies is in directory [E_assembly_evaluation](E_assembly_evaluation).

### Assembly evaluation

The choice of parameters for the final assembly was combination of several factors : *Continuity* (N50), *Completeness* (% of Ns, BUSCO) and *comparability* of sexual and asexual sister species.

The pipeline will calculate all the stats and pull them into table.

#### execution

Following four commands can be run in parallel. Calculate basic statistics:

```
make asm.stats
```

Search for single copy orthologs using BUSCO

```
make busco
```

Calculate number of unknown nucleotides

```
make counts.Ns
```

Once all the tasks above are done, pull the results into a table

TODO : pull them into table (I bet I got this already somewhere)
