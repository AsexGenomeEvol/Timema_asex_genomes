# Make stat files
include template.mk

## make asm.stats : will compute stats for all assemblies in ../data/*
.PHONY : asm.stats
asm.stats : ctg.stats scf.stats

## make ctg.stats : will compute stats for all contigs in ../data/*
.PHONY : ctg.stats
ctg.stats:
	$(MAKE) -s TYPE=ctg -f E_assembly_evaluation/individual_stats.mk

## make scf.stats : will compute stats for all scaffolds in ../data/*
.PHONY : scf.stats
scf.stats:
	$(MAKE) -s TYPE=scf -f E_assembly_evaluation/individual_stats.mk

CTG_LISTS = $(patsubst %, stats/assemblies/%_ctgs, $(SPECIES))
SCF_LISTS = $(patsubst %, stats/assemblies/%_scfs, $(SPECIES))

## make asm.tables : make all assembly tables (per species per ctg/scf)
.PHONY : asm.tables
asm.tables : ctg.tables scf.tables

## make ctg.tables : make tables of assemblies of contigs per species
.PHONY : ctg.tables
ctg.tables : $(CTG_LISTS)

## make scf.tables : make tables of assemblies of scaffolds per species
.PHONY : scf.tables
scf.tables : $(SCF_LISTS)

CMP_STAT_SCRIPT = E_assembly_evaluation/compute_stat_table.R

stats/assemblies/%_ctgs : $(CMP_STAT_SCRIPT)
	$(MAKE) -s SCRIPT=$< SP=$* TYPE=ctg -f E_assembly_evaluation/stats_lists.mk

stats/assemblies/%_scfs : $(CMP_STAT_SCRIPT)
	$(MAKE) -s SCRIPT=$< SP=$* TYPE=scf -f E_assembly_evaluation/stats_lists.mk

## make assemblies/ctgs_fulltable.tsv : make one big table of contig assemblies
## make assemblies/scfs_fulltable.tsv : make one big table of scaffold assemblies
assemblies/%_fulltable.tsv : E_assembly_evaluation/make_full_table.R stats/assemblies/*_%.tsv
	Rscript $< $*

print.parameter.effects.scf.% :
	Rscript E_assembly_evaluation/print_assembly_par_effects.R scf $*

print.parameter.effects.ctg.% :
	Rscript E_assembly_evaluation/print_assembly_par_effects.R ctg $*
