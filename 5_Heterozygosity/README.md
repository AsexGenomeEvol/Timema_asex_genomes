### UNDER CONSTRUCTION

1. [A_mapping_and_preprocessing](A_mapping_and_preprocessing)
2. [B_SNP_calling](B_SNP_calling)
3. [C_SV_calling](C_SV_calling)
4. [D_variant_analysis](D_variant_analysis)


Alright, all the SV calls take for ever detecting breakpoint between scaffolds as SVs (scaffolds are considered chromosomes in the SV world). To prevent this unintended behaviour I need to get rid of all reads mapping to edges of different scaffolds.
I was considering `BAMQL`, but it was too complicated to install on the cluster.
So, I will map `A_mapping_and_preprocessing/filter_splitreads.py`


TODO

```
5_Heterozygosity/B_SV_calling/coverages_table.R
5_Heterozygosity/B_SV_calling/prepare_coverage_tables.R
```
