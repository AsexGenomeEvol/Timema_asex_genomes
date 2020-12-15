# BATCH
include template.mk

ifeq ($(BATCH),1)
BUSCO_CHECKPOINTS := $(patsubst %, data/%/checkpoints/batch1_busco, $(SPECIES))
endif

ifeq ($(BATCH),2)
BUSCO_CHECKPOINTS := $(patsubst %, data/%/checkpoints/batch2_busco, $(SPECIES))
endif

ifeq ($(BATCH),3)
BUSCO_CHECKPOINTS := $(patsubst %, data/%/checkpoints/batch3_busco, $(SPECIES))
endif

.PHONY : busco
busco : $(BUSCO_CHECKPOINTS)

data/%/checkpoints/batch1_busco :
	$(eval KMER := $(shell Rscript scripts/R/get_optimal_k.R $* 1))
	cd data/$*/assembly/logs && \
		$$TROOT/E_assembly_evaluation/2_busco.sh $* \
			$$TROOT/data/$*/assembly/SOAP_k$(KMER)_nomse_nompe_fewdata/$*_SOAP_k$(KMER)_filt_GC.fa \
			SOAP_k$(KMER)_nomse_nompe_fewdata
	touch $@

data/%/checkpoints/batch2_busco :
	$(eval KMER := $(shell Rscript scripts/R/get_optimal_k.R $* 2))
	cd data/$*/assembly/logs && \
		$$TROOT/E_assembly_evaluation/2_busco.sh $* \
			$$TROOT/data/$*/assembly/SOAP_k$(KMER)_nomse_nompe/$*_SOAP_k$(KMER)_filt_GC.fa \
			SOAP_k$(KMER)_nomse_nompe
	touch $@

data/%/checkpoints/batch3_busco :
	$(eval KMER := $(shell Rscript scripts/R/get_optimal_k.R $*))
	cd data/$*/assembly/logs && \
		$$TROOT/E_assembly_evaluation/2_busco.sh $* \
			$$TROOT/data/$*/assembly/abyss_k$(KMER)_nomse_nompe_fewdata_BESST/Scaffolds_pass5_GC.fa \
			abyss_k$(KMER)_nomse_nompe_fewdata_BESST
	touch $@
