include template.mk

KMERGENIE := $(patsubst %, data/%/trimmed_reads/kmergenie, $(SPECIES))

.PHONY : kmergenie
kmergenie : $(KMERGENIE)

stats/reads/kmertable_filt.tsv : B_read_parsing/create_opt_kmer_table.R data/*/trimmed_reads/kmergenie
	Rscript $<

data/%/trimmed_reads/kmergenie :
	mkdir -p data/$*/trimmed_reads/logs;
	cd data/$*/trimmed_reads/logs ; \
		$$TROOT/B_read_parsing/kmergenie_lsf.sh $*

