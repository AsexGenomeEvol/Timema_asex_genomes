# TYPE - ctg (contigs) or scf (scaffolds)

ifeq ($(TYPE),ctg)
CTG_STATS = $(shell find data/$(SP) -type f -name '*_ctg.stats')
stats/assemblies/$(SP)_ctgs.tsv : $(SCRIPT) $(CTG_STATS)
	Rscript $^ $(SP)
else
SCF_STATS = $(shell find data/$(SP) -type f -name '*_scf.stats')
stats/assemblies/$(SP)_scfs.tsv : $(SCRIPT) $(SCF_STATS)
	Rscript $^ $(SP)
endif
