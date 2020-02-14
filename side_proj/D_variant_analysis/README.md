## SNV analysis

Extract genotypes out of heavy vcf files. Once we get the base stats we can come back to the original files to check other detailed info, but for now the genotypes only will do the job.

```bash
for sp  in $TIMEMAS; do
    grep "PASS" "$sp".SNP_filter.vcf | awk '{ print $1 "\t" $2 "\t" substr($10, 1, 3) "\t" substr($11, 1, 3) "\t" substr($12, 1, 3) "\t" substr($13, 1, 3) "\t" substr($14, 1, 3) }' > "$sp".SNP_filter_passed.tsv
done
```

This bash-prefiltered files I process with python script `sorting_variants.py`

```
python3 sorting_variants.py
```

to get three files: `*_heterozygous_SNP_filter_passed.tsv`, `*_homozygous_SNP_filter_passed.tsv`, and `*_trinalge_SNP_filter_passed.tsv`. The first two are lists of locations of SNPs found in at least two individuals that were found only in homozygous or only in heterozygous states. The triangle file are values that should be placed on the triangle plot (SFS decomposed by heterozygotes).

