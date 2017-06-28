/Home/kjaron/archive/timema_assembly/data/$(SP)/assembly/$(SP)_asm_b1.tar.gz :
	$(eval KMER := $(shell Rscript scripts/R/get_optimal_k.R $(SP) 1))
	echo archiving : data/$(SP)/assembly/SOAP_k$(KMER)_nomse_nompe_fewdata
	cd data/$(SP)/assembly && \
		tar cf - SOAP_k$(KMER)_nomse_nompe_fewdata | \
		pigz -9 -p 32 > $(SP)_asm_b1.tar.gz && \
		mv $(SP)_asm_b1.tar.gz $@  && \
		rm -r SOAP_k$(KMER)_nomse_nompe_fewdata
