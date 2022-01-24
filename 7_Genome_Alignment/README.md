
### Plotting variants on chromosomes

Panel B of Figure 3 can be regenerated using following script (source files are here in the repository):

```

scripts/plot_variants_on_chromosomes.R
```

### Regenerating source files for the plotting

The provided source files for the script above can be regenerated using several input files that would need to be computed beforehand, assuming

- `data/b3v08_anchoring_to_LGs/<SP>_scf_block_alignment.tsv`
- `data/SNP_calls/<SP>_reduced_filtered_variants.tsv` - generated in `5_Heterozygosity/SNP_calling`
- `data/manta_SV_calls/data/<SP>/SVs_filt_stringent_union.vcf` - generated in `6_Structural_Variants`

where `<SP>` is a placeholder for species codes (1_Tdi, 1_Tps, etc). If all these files are present, following script

```
Rscript scripts/make_variant_density_table.R
```

will generate a `tables/<SP>_variants_on_chromosomes_w1e+06.tsv` file for each species.

### Other (unpublished) plots

- `scripts/plot_SNP_density_variation.R` - plot of density variations of SNPs between species
- `colour_loci` - directory with expliration plots of the colour-encoding locus on lg8
