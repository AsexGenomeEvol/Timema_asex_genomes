BAM_FILES := $(shell find data -name "*.bam")
BAM_INDEXES := $(patsubst %, %.bai, $(BAM_FILES))

.PHONY : index.all
index.all : $(BAM_INDEXES)

%.bai : %
	$(eval DIR := $(shell dirname $*))
	$(eval BAM := $(shell basename $*))
	cd $(DIR) && $$TROOT/scripts/index_bam.sh $(BAM)
