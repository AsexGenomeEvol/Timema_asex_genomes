include template.mk

# VARIABLES

# names of checkpoint files (placeholders for submited jobs or multiple files)
TRIMMING_CHECKPOINTS := $(patsubst %, data/%/checkpoints/trimming, $(SPECIES))
VERIFIED_TRIMMING_CHECKPOINTS := $(patsubst %, trimming_verified, $(TRIMMING_CHECKPOINTS))
TRIMMED_READS_STATS := $(patsubst %, stats/reads/%_trimmed_reads.tsv, $(SPECIES))

# names of phantom targets
CLEAN_RAW_READS := $(patsubst %, clean.raw_reads.%, $(SPECIES))

# API

# make trimmed_reads : trimm all reads (require raw_reads)
.PHONY: trimmed_reads
trimmed_reads: $(TRIMMING_CHECKPOINTS)

# make clean.raw_reads : delete raw reads that are already trimmed
.PHONY : clean.raw_reads
clean.raw_reads : $(CLEAN_RAW_READS)

.PHONY : stats.trimmed_reads
stats.trimmed_reads : $(TRIMMED_READS_STATS)

# RULES

# TODO, add raw_reads checkpoint as dependency once I will have dl links
data/%/checkpoints/trimming :
	mkdir -p data/$*/trimmed_reads/logs data/$*/checkpoints
	cd data/$*/trimmed_reads/logs; \
	  $$TROOT/B_read_parsing/trim_pair_end_reads_lsf.sh $* is_350 1612 _run2; \
	  $$TROOT/B_read_parsing/trim_pair_end_reads_lsf.sh $* is_350 1602; \
	  $$TROOT/B_read_parsing/trim_pair_end_reads_lsf.sh $* is_550; \
	  $$TROOT/B_read_parsing/trim_pair_end_reads_lsf.sh $* is_700; \
	  $$TROOT/B_read_parsing/process_mate_pair_reads.sh $*
	touch $@


data/%/checkpoints/trimming_verified :
	bash B_read_parsing/verify_trimming.sh $* && touch $@

clean.raw_reads.% : data/%/checkpoints/trimming_verified
	rm -r data/$*/raw_reads/is_*

stats/reads/%_trimmed_reads.tsv : data/%/checkpoints/trimming_verified
	mkdir -p stats/reads
	find data/$*/trimmed_reads/ -name "*.fq.gz" -exec scripts/generic_genomics/fastq.gz2number_of_nt.sh {} \; > $@

stats/reads/coverage.table.tsv : B_read_parsing/theoretical_coverage.R $(TRIMMED_READS_STATS)
	Rscript $< $@
