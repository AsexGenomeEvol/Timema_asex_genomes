include template.mk

BATCH1_FILTERED := $(patsubst %, filt.batch1.%, $(SPECIES))
BATCH2_FILTERED := $(patsubst %, filt.batch2.%, $(SPECIES))
SPECIES_ARCHIVE.B1 := $(patsubst %, %.archive.b1, $(SPECIES))

## make filter.batch1 : filter contigs of batch 1 to greater than 250
.PHONY : filter.batch1
filter.batch1 : $(BATCH1_FILTERED)

## make filter.batch2 : filter contigs of batch 2 to greater than 250
.PHONY : filter.batch2
filter.batch2 : $(BATCH2_FILTERED)

## make archive.b1 : move b1 assemblies to archive
.PHONY : archive.b1
archive.b1 : $(SPECIES_ARCHIVE.B1)

filt.batch1.% :
	@echo $*
	$(MAKE) -s SP=$* BATCH=1 -f scripts/filter_asm.mk

filt.batch2.% :
	@echo $*
	$(MAKE) -s SP=$* BATCH=2 -f scripts/filter_asm.mk

%.archive.b1 :
	$(MAKE) /Home/kjaron/archive/timema_assembly/data/$*/assembly/$*_asm_b1.tar.gz -s SP=$* -f E_assembly_evaluation/filtering.mk

/Home/kjaron/archive/timema_assembly/data/$(SP)/assembly/$(SP)_asm_b1.tar.gz :
	$(eval KMER := $(shell Rscript scripts/R/get_optimal_k.R $(SP) 1))
	echo archiving : data/$(SP)/assembly/SOAP_k$(KMER)_nomse_nompe_fewdata
	cd data/$(SP)/assembly && \
		tar cf - SOAP_k$(KMER)_nomse_nompe_fewdata | \
		pigz -9 -p 32 > $(SP)_asm_b1.tar.gz && \
		mv $(SP)_asm_b1.tar.gz $@  && \
		rm -r SOAP_k$(KMER)_nomse_nompe_fewdata
