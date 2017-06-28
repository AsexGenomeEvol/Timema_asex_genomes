include template.mk
# VERSION - b3v04 for instance

CHECKPOINTS := $(foreach SP, $(SPECIES), data/$(SP)/checkpoints/ref_is350_to_$(VERSION))
REF_RG := "@RG\tID:ref\tLB:is350"

.PHONY : map.all
map.all : $(CHECKPOINTS)

# data/$*/mapping/ref_is350_to_$*_$(VERSION).bam
# data/1_Tdi/checkpoints/ref_is350_to_b3v04
data/%/checkpoints/ref_is350_to_$(VERSION) : 
	mkdir -p data/$*/mapping/logs
	cd data/$*/mapping/logs && \
		$$TROOT/scripts/map_ref_is350_RG_lsf.sh $* $(VERSION) $(REF_RG)
	touch $@
