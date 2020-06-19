## Analysis of SVs

Analysis of 6 individuals per species. I want to know heterozygosity patterns in individual species. Here we do:

 - basic population genetics on SVs - site frequency spectra, heterozygosity analysis
 - the density of SVs on scaffolds (addressing putative reproduction mode)
 - per type / per size analysis of the SVs in the genomes.

#### What do we have in the end

- individual delly / smove / manta SV calls (mostly if we want to go back and check something)
- merged calls of the three (`"$SP"_all_calls_merged.bcf`, is there support inside? Need to check)
- merged genotyping calls (`"$SP"_delly_genotyping_merged.bcf`), given the set of candidates

- merged paragraph genotyping (`data/genotyping/3_Tms_merged_calls_naive.vcf`), only for Tms for now. The merging procedure should be adjusted, therefore for now the file is called "naive".

### SV curation

This is sort of a post-analysis step to make sure what we do makes sense.

We will use [SV-plaudit](https://github.com/jbelyeu/SV-plaudit), a high throughput SV curation tool. It internally uses [samplot](https://github.com/ryanlayer/samplot) to plot individual SVs, I will try to use it first internally, then I will run the whole framework

```
samplot plot -n Tms_00 Tms_01 Tms_02 Tms_03 Tms_04 Tms_05 -b data/mapped_reseq_reads/Tms_00_to_b3v08_mapped_within_scfs.bam data/mapped_reseq_reads/Tms_01_to_b3v08_mapped_within_scfs.bam data/mapped_reseq_reads/Tms_02_to_b3v08_mapped_within_scfs.bam data/mapped_reseq_reads/Tms_03_to_b3v08_mapped_within_scfs.bam data/mapped_reseq_reads/Tms_04_to_b3v08_mapped_within_scfs.bam data/mapped_reseq_reads/Tms_05_to_b3v08_mapped_within_scfs.bam -o data/sandbox/3_Tms_b3v08_scaf000092_INV_80652_80811.png -s 80652 -e 80811 -c 3_Tms_b3v08_scaf000092 -a -t INV
```

this looks bad. Like really bad. Let's try to get them for all

```
qsub -o logs/ -e logs/ -cwd -N samplot -V -pe smp64 1 -b yes 'samplot vcf --vcf data/genotyping/3_Tms_merged_calls_naive_header.vcf -d figures/SVs -O png -b data/mapped_reseq_reads/Tms_00_to_b3v08_mapped_within_scfs.bam data/mapped_reseq_reads/Tms_01_to_b3v08_mapped_within_scfs.bam data/mapped_reseq_reads/Tms_02_to_b3v08_mapped_within_scfs.bam data/mapped_reseq_reads/Tms_03_to_b3v08_mapped_within_scfs.bam data/mapped_reseq_reads/Tms_04_to_b3v08_mapped_within_scfs.bam data/mapped_reseq_reads/Tms_05_to_b3v08_mapped_within_scfs.bam --sample_ids Tms_00 Tms_01 Tms_02 Tms_03 Tms_04 Tms_05 > F_structural_variation_analysis/samplot_commands.sh'
```

I should also add genes to [the plots](https://github.com/ryanlayer/samplot/tree/v1.0.1#gene-and-other-genomic-feature-annotations)

### basic population genetics

- SFS (site freq. spectra)
- SFS vs heterozygosity (figure out how to plot it, maybe monochromatic scale and a simple tile per type??)
- stats as fractions of "present in all", "heterozygous in all", "homozygous in all", "heterozygous"

Lot of the original plotting is done in very unsorted way in `F_structural_variation_analysis/playground.R` script using functions from `F_structural_variation_analysis/plot_SFS.R`.

The basic popgen of the final genotypes will now be done by

```
Rscript F_structural_variation_analysis/plot_genotyped_SV_SFS.R
```

There seems to be the last 12

```
# solid ones:
# 3_Tms_b3v08_scaf002312: 21725 - 23286; 00
# 3_Tms_b3v08_scaf000092: 80652 - 80811; 00 01 02 05
```

```
samtools view data/mapped_reseq_reads/Tms_00_to_b3v08_mapped_within_scfs.bam 3_Tms_b3v08_scaf002312 > data/mapped_reseq_reads/subselected_chromosomes/Tms_00_b3v08_scaf002312.bam
samtools view data/mapped_reseq_reads/Tms_01_to_b3v08_mapped_within_scfs.bam 3_Tms_b3v08_scaf002312 > data/mapped_reseq_reads/subselected_chromosomes/Tms_01_b3v08_scaf002312.bam
for i in 00 01 02 03 04 05; do
    samtools view data/mapped_reseq_reads/Tms_"$i"_to_b3v08_mapped_within_scfs.bam 3_Tms_b3v08_scaf000092 > data/mapped_reseq_reads/subselected_chromosomes/Tms_"$i"_b3v08_scaf000092.bam
done
```

### density of SVs

- check modality in sexual / asexuals (sexual are expected to segregate regardless of position in the genome, i.e. all scaffolds should be even) --> too few data, probably got to map them to the reference
- in asexuals we expect unimodal disr for homozygous SVs and bimodal for heterozygous SVs. Furthermore, we expect a negative relation of the two
- A statistical approach? Given a random placement, what is the expected distribution of observed density given the scf size, quantile? -> this could be useful

```
F_structural_variation_analysis/SV_density.R
```

Using all this, I can not detect any displacement of homozygous / heterozygous SV between asexuals. Supposedly this can be a problem of the combinations of too fragmented assemblies, no filtering applied to the called SVs and because SVs are rare and the signal was weak at the very beginning. I also calculated expecations of SV presence based on SV homo/heterozygosity, maybe I should figure out a sexual expececation jointly estimating expectations for both (that would more easily show differences between sex/asex)

### per Type / Size analysis

- For this I should probably separate homozygous and heterozygous SVs in asexuals (I will have 5 inds for homo, 6 for hetero).

```
F_structural_variation_analysis/SV_types_and_lengths.R
```

### SVs vs scaffold size where they were called

To figure:
 - How long scaffolds are meaningful to anchor to genome reference
 - What is the effect of assembly fragmentation on SV calling

```
F_structural_variation_analysis/SVs_vs_scaffold_length.R
```

### Sanity checks

I need to make sure that I understand how much these SV calls are wrong. I need to do a bunch of sanity checks, especially addressing coverage supporting individual SVs.

- Check coverage distribution of heterozygous alleles (per type?)
- Check coverage distribution of homozygous alleles (per type?)
- Table of coverages of individual libraries? (calculate from bam file?)