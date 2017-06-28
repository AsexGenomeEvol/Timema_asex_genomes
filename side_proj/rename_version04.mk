# SP

data/$(SP)/reference/$(SP)_b1v04.fa :
	$(eval KMER := $(shell Rscript scripts/R/get_optimal_k.R $(SP) 1))
	seqkit sort --by-length --reverse \
		data/$(SP)/assembly/SOAP_k$(KMER)_nomse_nompe_fewdata/$(SP)_SOAP_k$(KMER)_filt_GC.fa -o - \
		| scripts/rename_headers.awk -v sp=$(SP) version=b1v04 - > $@

data/$(SP)/reference/$(SP)_b3v04.fa :
	$(eval KMER := $(shell Rscript scripts/R/get_optimal_k.R $(SP) 3))
	seqkit sort --by-length --reverse \
		data/$(SP)/assembly/abyss_k$(KMER)_nomse_nompe_fewdata_BESST/Scaffolds_pass5_GC.fa -o - \
		| scripts/rename_headers.awk -v sp=$(SP) version=b3v04 - > $@
