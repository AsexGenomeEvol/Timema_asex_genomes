# Heterozygosity

Genome-wide nucleotide heterozygosity was estimated using genome profiling analysis of raw reads from the reference genomes using GenomeScope (v2). A second, SNP-based heterozygosity estimate was generated using re-sequenced individuals. We re-sequenced five individuals per species (and generated SNP and SV calls for all of them). However, 3 individuals of T. shepardi, 2 individuals of T. poppensis and one _T. tahoe_ individual did not pass quality control and were discarded from the downstream analyses.

- [Heterozygosity](#heterozygosity)
  * [Pipeline](#pipeline)
    + [Genome profiling](#genome-profiling)
    + [SNP calling](#snp-calling)
      - [execution](#execution)

## Pipeline

In this document we show how heterozygosity was estimated using genome profiling and SNPs.

### Genome profiling

**Genome_profiling** was used to generate genome wide estiamtes of heterozygosity. These analyses fail to distingust completelly homozygous genomes and genomes with very low heterozygosity levels. We used SNP calling to complement genome profiling analysis.

### SNP calling

**SNP calling** ([SNP_calling](SNP_calling)) was based on the GATK best practices pipeline. We used a conservative set of SNPs with quality scores ≥300, and supported by 15x coverage in at least one of the individuals. SNP heterozygosity was then estimated as the number of heterozygous SNPs divided by the number of callable sites in each genome. Due to stringent filtering criteria, our SNP based heterozygosity is an underestimation of genome-wide heterozygosity. For sexual species the GenomeScope estimate will be must closter to the biological reality.

#### execution

The SNP calling pipeline is completelly described in [SNP_calling](SNP_calling). The we get all the final output of SNPs calling in one directory: `data/SNP_calls/"$sp".SNP_filter.vcf.gz`, where `$sp` stands for the ten species codes. These output files will be used in two following sections.

Note all analysis of **SNP calls** are done together with **SV calls** ([6_Structural_Variants](../6_Structural_Variants)) and using **genome alignment** ([7_Genome_Alignment](../7_Genome_Alignment)) and are found in corresponding directories.
