
#### Pipeline steps for calling SNPs (and indels) in stick insects are described in :


* [**A - raw data formatting**](./A_raw_reads)

  Organising / merging / renaming fastq files of the five resequenced individuals in the ten species.

* [**B - read trimming**](./B_cleaned_reads)

  Trimming of the reads for quality and presence of adapters.

* [**C - mapping**](./C_mapping)

  Mapping trimmed reads on reference genomes.

* [**D - snp calling (preliminary round)**](./D_snp_calling_round0)

  First round of SNP (and indels) calling, this preliminary set of snps will only be used to mask positions during the following "base quality score recalibration" step. 

* [**E - base quality score recalibration (BQSR)**](./E_recalibration)

  Readjust base quality scores in bam files (obtained at the mapping step) by applying machine learning to empirically model systematic errors associated with each sequencing run.

* [**F - snp calling (final round)**](./F_snp_calling_round1)

  Second round of SNP (and indels) calling using updated bam files (with corrected base quality score).
