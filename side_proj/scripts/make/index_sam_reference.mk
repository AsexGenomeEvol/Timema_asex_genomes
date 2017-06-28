include template.mk
# SPECIES
# VERSION - b3v04 for instance

GENOMES := $(foreach SP, $(SPECIES), data/$(SP)/reference/$(SP)_$(VERSION).fa)
INDEXES := $(patsubst %, %.fai, $(GENOMES))

.PHONY : index.all
index.all : $(INDEXES)

%.fai : %
	$(eval DIR := $(shell dirname $*))
	$(eval GEN := $(shell basename $*))
	cd $(DIR) && $$VROOT/A_mapping_and_preprocessing/index_fa_sam.sh $(GEN)
