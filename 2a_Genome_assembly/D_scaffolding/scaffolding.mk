include template.mk

BATCH3_CONTIGS := $(shell find data/ -name "*-6.fa" | grep "nomse_nompe_fewdata")
BATCH3_FILTERED := $(patsubst %.fa, %_filt.fa, $(BATCH3_CONTIGS))
BATCH3_INDEX := $(patsubst %, data/%/checkpoints/index_batch3, $(SPECIES))
BATCH3_MAP := $(patsubst %, data/%/checkpoints/mapping_to_batch3, $(SPECIES))
BATCH3_BESST_CHECKPOINTS := $(patsubst %, data/%/checkpoints/batch3_BESST, $(SPECIES))
GAPFILLING_CHECKPOINTS1 := $(patsubst %, data/%/checkpoints/batch1_gc, $(SPECIES))
#GAPFILLING_CHECKPOINTS += $(patsubst %, data/%/checkpoints/batch2_gc, $(SPECIES))
GAPFILLING_CHECKPOINTS3 := $(patsubst %, data/%/checkpoints/batch3_gc, $(SPECIES))
BATCH3_REASM_INDEX := $(patsubst %, data/%/checkpoints/reasm_index_batch3, $(SPECIES))
BATCH3_REASM_MAP := $(patsubst %, data/%/checkpoints/reasm_mapping_to_batch3, $(SPECIES))
BATCH3_REASM_BESST_CHECKPOINTS := $(patsubst %, data/%/checkpoints/reasm_batch3_BESST, $(SPECIES))
BATCH3_REASM_GAPFILLING := $(patsubst %, data/%/checkpoints/reasm_batch3_gc, $(SPECIES))

## make filter.batch3 : filter contigs of batch 3 to greater than 250
.PHONY : filter.batch3
filter.batch3 : $(BATCH3_FILTERED)

## make indexed.batch3 : index filtered contigs of batch 3
.PHONY : indexed.batch3
indexed.batch3 : $(BATCH3_INDEX)

.PHONY : indexed.reasm.batch3
indexed.reasm.batch3 : $(BATCH3_REASM_INDEX)

## make mapping.batch3 : map pe reads to contig assembly of batch3
.PHONY : mapping.batch3
mapping.batch3 : $(BATCH3_MAP)

.PHONY : mapping.reasm.batch3
mapping.reasm.batch3 : $(BATCH3_REASM_MAP)

## make scaffold.batch3 : scaffold abyss contigs with BESST
.PHONY : scaffold.batch3
make scaffold.batch3 : $(BATCH3_BESST_CHECKPOINTS)

.PHONY : scaffold.reasm.batch3
scaffold.reasm.batch3 : $(BATCH3_REASM_BESST_CHECKPOINTS)

## make gapfilling
.PHONY : gapfilling
gapfilling : $(GAPFILLING_CHECKPOINTS1) $(GAPFILLING_CHECKPOINTS3)

.PHONY : gapfilling.batch1
gapfilling.batch1 : $(GAPFILLING_CHECKPOINTS1)

.PHONY : gapfilling.batch3
gapfilling.batch3 : $(GAPFILLING_CHECKPOINTS3)

.PHONY : gapfilling.reasm.batch3
gapfilling.reasm.batch3 : $(BATCH3_REASM_GAPFILLING)

%_filt.fa : %.fa
	D_scaffolding/contig_filtering.sh $< $@

data/%/checkpoints/index_batch3 :
	$(eval KMER := $(shell Rscript scripts/R/get_optimal_k.R $*))
	$(eval DIR := data/$*/assembly/abyss_k$(KMER)_nomse_nompe_fewdata)
	$(eval REF := $*_abyss_k$(KMER)-6_filt.fa)
	cd $(DIR) && $$TROOT/scripts/index_fa_bwa.sh $(REF)
	touch $@

data/%/checkpoints/mapping_to_batch3 :
	$(eval KMER := $(shell Rscript scripts/R/get_optimal_k.R $*))
	$(eval DIR := data/$*/assembly/abyss_k$(KMER)_nomse_nompe_fewdata/mapping)
	$(eval REF := ../$*_abyss_k$(KMER)-6_filt.fa)
	$(eval RDIR := ../../../trimmed_reads)
	mkdir -p $(DIR)
	cd $(DIR); \
	  $$TROOT/D_scaffolding/map_reads_for_BESST.sh $*_350 $(REF) $(RDIR)/$*_R1t_is_350.fq.gz $(RDIR)/$*_R2t_is_350.fq.gz ; \
	  $$TROOT/D_scaffolding/map_reads_for_BESST.sh $*_550 $(REF) $(RDIR)/$*_R1t_is_550.fq.gz $(RDIR)/$*_R2t_is_550.fq.gz ; \
	  $$TROOT/D_scaffolding/map_reads_for_BESST.sh $*_700 $(REF) $(RDIR)/$*_R1t_is_700.fq.gz $(RDIR)/$*_R2t_is_700.fq.gz ; \
	  $$TROOT/D_scaffolding/map_reads_for_BESST.sh $*_3000 $(REF) $(RDIR)/mp_nxtrim_FR/$*_R1t_is_3000.fq.gz $(RDIR)/mp_nxtrim_FR/$*_R2t_is_3000.fq.gz ; \
	  $$TROOT/D_scaffolding/map_reads_for_BESST.sh $*_5000 $(REF) $(RDIR)/mp_nxtrim_FR/$*_R1t_is_5000.fq.gz $(RDIR)/mp_nxtrim_FR/$*_R2t_is_5000.fq.gz ;
	touch $@

data/%/checkpoints/batch3_BESST :
	$(eval KMER := $(shell Rscript scripts/R/get_optimal_k.R $*))
	mkdir -p data/$*/assembly/logs;
	cd data/$*/assembly/logs; \
	  $$TROOT/D_scaffolding/1_BESST_lsf.sh $* abyss_k$(KMER)_nomse_nompe_fewdata $*_abyss_k$(KMER)-6_filt.fa
	touch $@

data/%/checkpoints/batch1_gc :
	$(eval KMER := $(shell Rscript scripts/R/get_optimal_k.R $* 1))
	cd data/$*/assembly/logs && \
		$$TROOT/D_scaffolding/2_gap_filling_lsf.sh $* \
			$$TROOT/data/$*/assembly/SOAP_k$(KMER)_nomse_nompe_fewdata/$*_SOAP_k$(KMER)_filt.fa \
			SOAP_k$(KMER)_nomse_nompe_fewdata
	touch $@

data/%/checkpoints/batch2_gc :
	$(eval KMER := $(shell Rscript scripts/R/get_optimal_k.R $* 2))
	cd data/$*/assembly/logs && \
		$$TROOT/D_scaffolding/2_gap_filling_lsf.sh $* \
			$$TROOT/data/$*/assembly/SOAP_k$(KMER)_nomse_nompe/$*_SOAP_k$(KMER)_filt.fa \
			SOAP_k$(KMER)_nomse_nompe
	touch $@

data/%/checkpoints/batch3_gc :
	$(eval KMER := $(shell Rscript scripts/R/get_optimal_k.R $*))
	cd data/$*/assembly/logs && \
        $$TROOT/D_scaffolding/2_gap_filling_lsf.sh $* \
            $$TROOT/data/$*/assembly/abyss_k$(KMER)_nomse_nompe_fewdata_BESST/BESST_output/pass5/Scaffolds_pass5.fa \
            abyss_k$(KMER)_nomse_nompe_fewdata_BESST
	touch $@

data/%/checkpoints/reasm_index_batch3 :
	$(eval KMER := $(shell Rscript scripts/R/get_optimal_k.R $* kmergenie.filt))
	$(eval DIR := data/$*/assembly/abyss_k$(KMER)_nomse_nompe_fewdata_filtered)
	$(eval REF := $*_abyss_k$(KMER)-6_filt.fa)
	cd $(DIR) && $$TROOT/scripts/index_fa_bwa.sh $(REF)
	touch $@

data/%/checkpoints/reasm_mapping_to_batch3 :
	$(eval KMER := $(shell Rscript scripts/R/get_optimal_k.R $* kmergenie.filt))
	$(eval DIR := data/$*/assembly/abyss_k$(KMER)_nomse_nompe_fewdata_filtered/mapping)
	$(eval REF := ../$*_abyss_k$(KMER)-6_filt.fa)
	$(eval RDIR := ../../../trimmed_reads)
	mkdir -p $(DIR)
	cd $(DIR); \
	  $$TROOT/D_scaffolding/map_reads_for_BESST.sh $*_350 $(REF) $(RDIR)/$*_R{1,2}t_is_350.no_contaminant.fastq.gz ; \
	  $$TROOT/D_scaffolding/map_reads_for_BESST.sh $*_550 $(REF) $(RDIR)/$*_R{1,2}t_is_550.no_contaminant.fastq.gz ; \
	  $$TROOT/D_scaffolding/map_reads_for_BESST.sh $*_700 $(REF) $(RDIR)/$*_R{1,2}t_is_700.no_contaminant.fastq.gz ; \
	  $$TROOT/D_scaffolding/map_reads_for_BESST.sh $*_3000 $(REF) $(RDIR)/mp_nxtrim_FR/$*_R{1,2}t_is_3000.no_contaminant.fastq.gz ; \
	  $$TROOT/D_scaffolding/map_reads_for_BESST.sh $*_5000 $(REF) $(RDIR)/mp_nxtrim_FR/$*_R{1,2}t_is_5000.no_contaminant.fastq.gz ;
	touch $@

data/%/checkpoints/reasm_batch3_BESST :
	$(eval KMER := $(shell Rscript scripts/R/get_optimal_k.R $* kmergenie.filt))
	mkdir -p data/$*/assembly/logs;
	cd data/$*/assembly/logs; \
	  $$TROOT/D_scaffolding/1_BESST_lsf.sh $* abyss_k$(KMER)_nomse_nompe_fewdata_filtered $*_abyss_k$(KMER)-6_filt.fa
	touch $@

data/%/checkpoints/reasm_batch3_gc :
	$(eval KMER := $(shell Rscript scripts/R/get_optimal_k.R $* kmergenie.filt))
	cd data/$*/assembly/logs && \
        $$TROOT/D_scaffolding/2_gap_filling_lsf.sh $* \
            $$TROOT/data/$*/assembly/abyss_k$(KMER)_nomse_nompe_fewdata_filtered_BESST/BESST_output/pass5/Scaffolds_pass5.fa \
            abyss_k$(KMER)_nomse_nompe_fewdata_BESST --filtered
	touch $@
