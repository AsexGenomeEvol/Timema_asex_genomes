include template.mk

SOAP_BATCH1k41_CHECKPOINTS := $(patsubst %, data/%/checkpoints/asm_batch1_k41, $(SPECIES))
SOAP_BATCH1k43_CHECKPOINTS := $(patsubst %, data/%/checkpoints/asm_batch1_k43, $(SPECIES))
SOAP_BATCH1k45_CHECKPOINTS := $(patsubst %, data/%/checkpoints/asm_batch1_k45, $(SPECIES))
BATCH1_OPT := $(patsubst %, data/%/checkpoints/asm_batch1_optimal, $(SPECIES))
BATCH2_OPT := $(patsubst %, data/%/checkpoints/asm_batch2_optimal, $(SPECIES))
BATCH3_CHECKPOINTS := $(patsubst %, data/%/checkpoints/asm_batch3, $(SPECIES))
BATCH3_REASSEMBLED_CHECKPOINTS := $(patsubst %, data/%/checkpoints/reasm_batch3, $(SPECIES))
CLEAN_ASM := $(patsubst %, clean.asm.%, $(SPECIES))

## make assemble.batch1 : submit a batch of assemblies (SOAP, kmers 41, 43, 45, --nomse --nompe --fewdata)
.PHONY : assemble.batch1
assemble.batch1 : $(SOAP_BATCH1k41_CHECKPOINTS) $(SOAP_BATCH1k43_CHECKPOINTS) $(SOAP_BATCH1k45_CHECKPOINTS)

## make assemble.batch1.next_iter : do next iteration in otpimal asm search (SOAP, --nomse --nompe --fewdata)
.PHONY : assemble.batch1.next_iter
assemble.batch1.next_iter : $(BATCH1_OPT)

## make assemble.batch2.next_iter : do an iteration in otpimal asm search (SOAP, --nomse --nompe)
.PHONY : assemble.batch2.next_iter
assemble.batch2.next_iter : $(BATCH2_OPT)

## make assemble.batch3 : submit a batch of assemblies (abyss, opt kmers, --nomse --nompe --fewdata)
.PHONY : assemble.batch3
assemble.batch3 : $(BATCH3_CHECKPOINTS)

## make clean.asm 2> stats/assemblies/current_optimal.tsv : clean all intermediate files from assmeblies that are suboptimal
.PHONY : clean.asm
clean.asm : $(CLEAN_ASM)

.PHONY : reassemble.batch3
reassemble.batch3 : $(BATCH3_REASSEMBLED_CHECKPOINTS)

data/%/checkpoints/asm_batch1_k41 :
	mkdir -p data/$*/assembly/logs;
	cd data/$*/assembly/logs; \
	$$TROOT/C_contig_assembly/2_SOAP.sh $* 41 --nomse --nompe --fewdata
	touch $@

data/%/checkpoints/asm_batch1_k43 :
	mkdir -p data/$*/assembly/logs;
	cd data/$*/assembly/logs; \
	$$TROOT/C_contig_assembly/2_SOAP.sh $* 43 --nomse --nompe --fewdata
	touch $@

data/%/checkpoints/asm_batch1_k45 :
	mkdir -p data/$*/assembly/logs;
	cd data/$*/assembly/logs; \
	$$TROOT/C_contig_assembly/2_SOAP.sh $* 45 --nomse --nompe --fewdata
	touch $@

data/%/checkpoints/asm_batch1_optimal :
	if [[ $$(Rscript scripts/R/submit_next_asm_batch.R $* 1) -eq 0 ]]; then touch $@; fi

data/%/checkpoints/asm_batch2_optimal :
	if [[ $$(Rscript scripts/R/submit_next_asm_batch.R $* 2) -eq 0 ]]; then touch $@; fi

data/%/checkpoints/asm_batch3 :
	$(eval KMER := $(shell Rscript scripts/R/get_optimal_k.R $*))
	mkdir -p data/$*/assembly/logs;
	cd data/$*/assembly/logs; \
	$$TROOT/C_contig_assembly/1_abyss.sh $* $(KMER) --nomse --nompe --fewdata
	touch $@

clean.asm.% :
	@echo $*
	cd stats; $(MAKE) assemblies/$*_ctgs.tsv assemblies/$*_scfs.tsv
	$(MAKE) -s SP=$* -f scripts/clean_asm.mk

data/%/checkpoints/reasm_batch3 :
	$(eval KMER := $(shell Rscript scripts/R/get_optimal_k.R $* kmergenie.filt))
	mkdir -p data/$*/assembly/logs;
	cd data/$*/assembly/logs; \
	$$TROOT/C_contig_assembly/1_abyss.sh $* $(KMER) --nomse --nompe --fewdata --filtered
	touch $@
