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

to get three files: `*_heterozygous_SNP_filter_passed.tsv`, `*_homozygous_SNP_filter_passed.tsv`, and `*_trinalge_SNP_filter_passed.tsv`. The first two are lists of locations of SNPs found in at least two individuals that were found only in homozygous or only in heterozygous states. The triangle file are values that should be placed on the triangle plot (SFS decomposed by heterozygotes).

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

