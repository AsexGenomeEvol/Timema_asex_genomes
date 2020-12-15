# BATCH
include template.mk

ifeq ($(BATCH),1)
QUAST_CHECKPOINTS := $(patsubst %, data/%/checkpoints/batch1_quast, $(SPECIES))
endif

ifeq ($(BATCH),2)
QUAST_CHECKPOINTS := $(patsubst %, data/%/checkpoints/batch2_quast, $(SPECIES))
endif

ifeq ($(BATCH),3)
QUAST_CHECKPOINTS := $(patsubst %, data/%/checkpoints/batch3_quast, $(SPECIES))
endif

.PHONY : quast
quast : $(QUAST_CHECKPOINTS)

data/%/checkpoints/batch1_quast :
	$(eval KMER := $(shell Rscript scripts/R/get_optimal_k.R $* 1))
	cd data/$*/assembly/logs && \
		$$TROOT/E_assembly_evaluation/1_quast.sh $* \
			$$TROOT/data/$*/assembly/SOAP_k$(KMER)_nomse_nompe_fewdata/$*_SOAP_k$(KMER)_filt_GC.fa \
			SOAP_k$(KMER)_nomse_nompe_fewdata
	touch $@

data/%/checkpoints/batch2_quast :
	$(eval KMER := $(shell Rscript scripts/R/get_optimal_k.R $* 2))
	cd data/$*/assembly/logs && \
		$$TROOT/E_assembly_evaluation/1_quast.sh $* \
			$$TROOT/data/$*/assembly/SOAP_k$(KMER)_nomse_nompe/$*_SOAP_k$(KMER)_filt_GC.fa \
			SOAP_k$(KMER)_nomse_nompe
	touch $@

data/%/checkpoints/batch3_quast :
	$(eval KMER := $(shell Rscript scripts/R/get_optimal_k.R $*))
	cd data/$*/assembly/logs && \
		$$TROOT/E_assembly_evaluation/1_quast.sh $* \
			$$TROOT/data/$*/assembly/abyss_k$(KMER)_nomse_nompe_fewdata_BESST/Scaffolds_pass5_GC.fa \
			abyss_k$(KMER)_nomse_nompe_fewdata_BESST
	touch $@
