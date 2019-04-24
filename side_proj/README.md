### Analysis of resequencing data

this repo will contain only analysis of resequencing data.

The analysis of variants of reference genomes were moved back to repository `timema_assembly`


#### Inference of variants

- SNP calling
- SV calling


#### Downstream analysis of variants

- ?
- ??
- ???

#### Dependencies

- [B_variant_calling](B_variant_calling) -
- [E_structural_variation_calling](E_structural_variation_calling) - `Delly` (0.7.6), `manta` (1.5.0)


-> you need to have genome reference that were build in the repository `timema_assembly`:

```
for sp in $TIMEMAS; do
    ln -s /scratch/beegfs/monthly/kjaron/timema_assembly/data/$sp/reference /scratch/beegfs/monthly/kjaron/variant_analysis/data/$sp/reference;
done
```

-> trimmed reads from Montpellier (check [A_mapping_and_preprocessing/README.md] for sorting them out)