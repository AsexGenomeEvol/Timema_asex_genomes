# Heterozygosity

Genome-wide nucleotide heterozygosity was estimated using genome profiling analysis of raw reads from the reference genomes using GenomeScope (v2). A second, SNP-based heterozygosity estimate was generated using re-sequenced individuals. We re-sequenced five individuals per species (and generated SNP and SV calls for all of them). However, 3 individuals of T. shepardi, 2 individuals of T. poppensis and one T. tahoe individual did not pass quality control and were discarded from the downstream analyses.



## Pipeline

TODO: table of content

### Genome profiling

**Genome_profiling** was used to generate genome wide estiamtes of heterozygosity. These analyses fail to distingust completelly homozygous genomes and genomes with very low heterozygosity levels. We used SNP calling to complement genome profiling analysis.

#### execution

### SNP calling

**SNP calling** ([A_SNP_calling](A_SNP_calling)) was based on the GATK best practices pipeline. We used a conservative set of SNPs with quality scores ≥300, and supported by 15x coverage in at least one of the individuals. SNP heterozygosity was then estimated as the number of heterozygous SNPs divided by the number of callable sites in each genome. Due to stringent filtering criteria, our SNP based heterozygosity is an underestimation of genome-wide heterozygosity. For sexual species the GenomeScope estimate will be must closter to the biological reality.

#### execution

The pipeline is completelly described in [A_SNP_calling](A_SNP_calling). Some addenitional filtering steps were deployed latter in [TODO](TODO) section.

### Visualization

TODO

#### Input data

**SNPs**

`data/SNP_calls/"$sp".SNP_filter.vcf.gz`
