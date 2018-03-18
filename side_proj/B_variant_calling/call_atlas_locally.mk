# MITE expected input

ifndef MITE
MITE := os5
endif

LIB := 350
WINDOW := 10000

GENOME := $(shell find $(MITE) -name "*_b1v03.fasta")
BAMS := $(wildcard $(MITE)/*$(LIB).sorted.bam)

MERGED_BAM := $(MITE)/$(MITE).trim.all.$(LIB).sorted.recal.bam
THETA_CALL := $(MITE)/$(MITE).trim.all.$(LIB).sorted.recal_theta_estimates.txt

RECAL_FILES := $(patsubst %.bam, %_recal.txt, $(BAMS))
RBAMS := $(patsubst %.bam, %.recal.bam, $(BAMS))

### calling targets

.PHONY : all
all : $(THETA_CALL)

### recepies

%.bam.bai : %.bam
        samtools index $<

.PRECIOUS : %_recal.txt
%_recal.txt : %.bam %.bam.bai $(GENOME)
        lacer.pl -bam $< -reference $(GENOME) \
                -outfile $@ -stopbases 3000000 -threads 4

%.recal.bam : %_recal.txt %.bam

.PRECIOUS : $(MERGED_BAM)
$(MERGED_BAM) : $(RBAMS)
        samtools merge $@ $^

$(THETA_CALL) : $(MERGED_BAM)
        atlas task=estimateTheta \
                bam=$< \
                window=$(WINDOW) \
                limitChr=$(shell python3 create_BED_withmin_window_size.py $(GENOME) $(WINDOW)) \
                suppressWarnings verbose \
                1> $(MITE)/theta.log
