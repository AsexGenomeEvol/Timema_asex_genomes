## SNV analysis

Extract genotypes out of heavy vcf files. Once we get the base stats we can come back to the original files to check other detailed info, but for now the genotypes only will do the job.

```bash
for sp in $TIMEMAS; do
    zcat data/SNP_calls/"$sp".SNP_filter.vcf.gz | grep -v "^#" | awk '{ if ($7 == "PASS"){ print $0 } }' |  cut -f 1,2,6,10,11,12,13,14 > data/SNP_calls/"$sp".SNP_filter_passed.tsv
done
```

This bash-prefiltered files I process with python script `sorting_variants.py`

```
python3 D_variant_analysis/sorting_variants.py
```

Produces `data/SNP_calls/<sp>_reduced_variants.tsv` which is just a reduced vcf files to release a bit of pain when loading them in R. The file contains variants, their corresponding scf and position, variant quality score, all 5 genotypes and corresponding depths supporting the genotypes.

### Analysing coverage of heterozygous SNPs

Variants have following INFO tags

```
GT:AD:DP:FT:GQ:PL

GT - Genotype 0/0 0/1 1/1
AD - Allelic depths for the ref and alt alleles in the order listed (separated by ,)
DP - Approximate read depth (reads with MQ=255 or with bad mates are filtered
FT - Genotype-level filter
GQ - Genotype Quality
PL - Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification
```

The script above `sorting_variants.py` produces also `_homo_depths.tsv` and `_hetero_depths.tsv` files with depths of variants that are found only in homozygous or only in heterozygous states.

In the script `process_variant_depths.R`, we plot a distribution of covrages per individual and a distribution of qualities of variants. In both cases we can see mulimodality. The most meaningful will be to filter all variants with at least one individual with coverage > 15 and variant quality score > 300.

### Filtering

I added this filtering criteria (at least one individual with coverage > 15 and variant quality score > 300) to the preprocessing python script, and now

```
python3 D_variant_analysis/sorting_variants.py
```

also generates `data/SNP_calls/<sp>_reduced_filtered_variants.tsv` files.

We could get some filtering stats I guess by `wc -l` ing the `_reduced` and `_reduced_filtered` files or by adding a few lines to the `sorting_variants` script.

### Plotting variants on chromosomes

I am resolving the plotting in the following script

```
Rscript D_variant_analysis/plot_SNPs_on_chromosomes.R
```

The files I need are called SNPs - `data/SNP_calls/<sp>_reduced_filtered_variants.tsv`, block alignment `data/b3v08_anchoring_to_LGs/<sp>_scf_block_alignment.tsv`, and genome index `data/<sp>/reference/<sp>_b3v08_scf.lengths`

