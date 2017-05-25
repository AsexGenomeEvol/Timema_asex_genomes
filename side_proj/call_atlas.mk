BAM_FILES := $(shell find data/*/mapping -name "*bam")
VARCALLS := $(patsubst %, %.varcall, $(BAM_FILES))

.PHONY : call.all
call.all : $(VARCALLS)

%.varcall :
	$(eval SP := $(shell echo $* | cut -f 2 -d /))
	mkdir -p data/$(SP)/variant_calls/logs
	$(eval REF := $(shell echo $* | cut -f 5 -d "_" | cut -f 1 -d '.'))
	$(eval IND := $(shell basename $* _to_$(REF).bam))
	$(eval MLE := data/$(SP)/variant_calls/$(IND)/atlas/$(IND)_to_$(REF)_MLEGenotypes.vcf.gz)
	$(eval THETA := data/$(SP)/variant_calls/$(IND)/atlas/$(IND)_to_$(REF)_theta_estimates.txt)
	$(MAKE) $(MLE) -s SP=$(SP) REF=$(REF) IND=$(IND) -f N_variant_calling/call_atlas.mk
	$(MAKE) $(THETA) -s SP=$(SP) REF=$(REF) IND=$(IND) -f N_variant_calling/call_atlas.mk

data/$(SP)/variant_calls/$(IND)/atlas/$(IND)_to_$(REF)_MLEGenotypes.vcf.gz :
	cd data/$(SP)/variant_calls/logs && \
		$$TROOT/N_variant_calling/atlas_callMLE.sh $(SP) $(REF) $(IND)

data/$(SP)/variant_calls/$(IND)/atlas/$(IND)_to_$(REF)_theta_estimates.txt :
	cd data/$(SP)/variant_calls/logs && \
		$$TROOT/N_variant_calling/atlas_estTheta.sh $(SP) $(REF) $(IND)
