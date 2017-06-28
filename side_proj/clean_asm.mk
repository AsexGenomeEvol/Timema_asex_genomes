DIR_TO_CLEAN := $(shell Rscript scripts/R/get_suboptimal_asm.R $(SP))
DIR_TARGETS := $(patsubst %, clean.asm.%, $(DIR_TO_CLEAN))

.PHONY : clean.asm
clean.asm : $(DIR_TARGETS)

clean.asm.% : data/$(SP)/assembly/%
	echo Cleaning $*;
	# abyss cleaning
	find $< -name "*-[123457].fa" -exec rm {} \;
	find $< -name "*-[123457].dot" -exec rm {} \;
	# SOAP cleaning
	find $< -name "*.updated.edge" -exec rm {} \;
	find $< -name "*readInGap.gz" -exec rm {} \;
	find $< -name "*readOnContig.gz" -exec rm {} \;
	find $< -name "*.links" -exec rm {} \;
	find $< -name "*.preArc" -exec rm {} \;
	find $< -name "*.edge.gz" -exec rm {} \;
	find $< -name "*.vertex" -exec rm {} \;
	find $< -name "*.path" -exec rm {} \;
	find $< -name "*.Arc" -exec rm {} \;
	find $< -name "*.newContigIndex" -exec rm {} \;
	find $< -name "*.markOnEdge" -exec rm {} \;
	find $< -name "*.ContigIndex" -exec rm {} \;
	find $< -name "*.scaf" -exec rm {} \;
