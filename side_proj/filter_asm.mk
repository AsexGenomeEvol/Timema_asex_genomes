# SP
# BATCH
KMER := $(shell Rscript scripts/R/get_optimal_k.R $(SP) $(BATCH))

ifeq ($(BATCH),1)
ASM = data/$(SP)/assembly/SOAP_k$(KMER)_nomse_nompe_fewdata/$(SP)_SOAP_k$(KMER).scafSeq
FILT = data/$(SP)/assembly/SOAP_k$(KMER)_nomse_nompe_fewdata/$(SP)_SOAP_k$(KMER)_filt.fa
else
ASM = data/$(SP)/assembly/SOAP_k$(KMER)_nomse_nompe/$(SP)_SOAP_k$(KMER).scafSeq
FILT = data/$(SP)/assembly/SOAP_k$(KMER)_nomse_nompe/$(SP)_SOAP_k$(KMER)_filt.fa
endif

$(FILT) : $(ASM)
	D_scaffolding/contig_filtering.sh $< $@
