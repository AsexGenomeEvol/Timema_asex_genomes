### SV calling

Complementing nucleotide-based measures we searched for **Structural variants**, to investigate if any heterozygosity is carried by SVs. We used Manta (v1.5.0), a diploid-aware pipeline for structural variant (SV) calling, in the same set of re-sequenced individuals used for SNP heterozygosity estimates.

We found a high frequency of heterozygous SVs with approximately twice the expected coverage (SM Figure 7), which likely represent false positives. To reduce the number of false positives, we filtered very short SVs (30 bases or less) and  kept only variant calls that had either split read or paired-end read support within the expected coverage range, where the coverage range was defined individually for each sample by manual inspection of coverage distributions. The filtered SV calls were subsequently merged into population SV calls using SURVIVOR (v1.0.2). The merging criteria were: SV calls of the same type on the same strand with breakpoints distances shorter than 100 bp.  



 -  [scripts](scripts)

 -  [C_variant_analysis](C_variant_analysis)


Alright, all the SV calls take for ever detecting breakpoint between scaffolds as SVs (scaffolds are considered chromosomes in the SV world). To prevent this unintended behaviour I need to get rid of all reads mapping to edges of different scaffolds.
I was considering `BAMQL`, but it was too complicated to install on the cluster.
So, I will map `A_mapping_and_preprocessing/filter_splitreads.py`


TODO

```
5_Heterozygosity/scripts/coverages_table.R
5_Heterozygosity/scripts/prepare_coverage_tables.R
```
