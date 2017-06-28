BAM_FILES := $(shell find data/*/mapping -name "*bam")
VARCALLS := $(patsubst %, %.varcall, $(BAM_FILES))

.PHONY : call.all
call.all : $(VARCALLS)

%.varcall :
	$(eval SP := $(shell echo $* | cut -f 2 -d /))
	mkdir -p data/$(SP)/variant_calls/logs
	$(eval REF := $(shell echo $* | cut -f 5 -d "_" | cut -f 1 -d '.'))
	$(eval IND := $(shell basename $* _to_$(REF).bam))
	$(eval MANTA := data/$(SP)/variant_calls/$(IND)/manta)
	$(eval DELLY := data/$(SP)/variant_calls/$(IND)/delly)
	$(MAKE) $(MANTA) -s SP=$(SP) REF=$(REF) IND=$(IND) \
		-f M_structural_variations/call_structural_variations.mk
	$(MAKE) $(DELLY) -s SP=$(SP) REF=$(REF) IND=$(IND) \
		-f M_structural_variations/call_structural_variations.mk

data/$(SP)/variant_calls/$(IND)/manta :
	cd data/$(SP)/variant_calls/logs &&
		$$TROOT/M_structural_variations/manta.sh $(SP) $(REF) $(IND)
	mkdir -p $@

data/$(SP)/variant_calls/$(IND)/delly :
	cd data/$(SP)/variant_calls/logs &&
		$$TROOT/M_structural_variations/delly.sh $(SP) $(REF) $(IND)
	mkdir -p $@
