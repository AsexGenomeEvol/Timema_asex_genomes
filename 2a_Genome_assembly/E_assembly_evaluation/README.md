## Testing and details

Evaluation of quality of assembly was based mainly on continuity (NG50), completeness (BUSCO score) and finally difference between sexual a asexual assembly (that we tried to minimize).

### Computing stats of assemblies

The stats will be calculated on all the assembles using (finds all assemblies and calculate stats)

```
make asm.stats
```

than you can create twenty tables of assemblies (two per species, one of contigs one of scaffilds)

```
make asm.tables
```

to create overall summary of contigs/scaffolds of all the assemblies, run

```
make stats/assemblies/ctgs_fulltable.tsv -f assembly_benchmarking_stats.mk
make stats/assemblies/scfs_fulltable.tsv -f assembly_benchmarking_stats.mk
```

### Reapr

I had several issues running reapr, therefore it will not be used for now (maybe latter on together with corrections using  structural variations).

### Quast

I also calculated QC with quast (number of ORFs, more detailed statistics). But in fact I did not use the output for anything, at the time there was no reasonable reference to compare it with

```
make quast
```
