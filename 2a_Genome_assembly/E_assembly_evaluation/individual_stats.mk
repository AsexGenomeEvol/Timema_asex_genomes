## Construct that ensures that files are searched only in case that stats supposed to be updated

ASM_PATH = data/*/assembly

ifeq ($(TYPE),ctg)
# find all ctg files
ABYSS_CTG = $(shell find $(ASM_PATH) -name '*-6.fa')
SOAP_CTG = $(shell find $(ASM_PATH) -name '*.contig')
# make one list of .stats files
CTG_STATS_TOCOMP = $(patsubst %.fa, %_ctg.stats, $(ABYSS_CTG))
CTG_STATS_TOCOMP += $(patsubst %.contig, %_ctg.stats, $(SOAP_CTG))
# make all ctg stats
.PHONY : ctg.stats
ctg.stats: $(CTG_STATS_TOCOMP)
else
# find all scaffold files
ABYSS_SCF = $(shell find $(ASM_PATH) -name '*-8.fa')
SOAP_SCF = $(shell find $(ASM_PATH) -name '*.scafSeq')
BESST_SCF = $(shell find $(ASM_PATH) -name 'Scaffolds_pass5.fa')

# make one stats list out of them
SCF_STATS_TOCOMP = $(patsubst %.fa, %_scf.stats, $(ABYSS_SCF))
SCF_STATS_TOCOMP += $(patsubst %.scafSeq, %_scf.stats, $(SOAP_SCF))
SCF_STATS_TOCOMP += $(patsubst %.fa, %_scf.stats, $(BESST_SCF))
# make all scf stats
.PHONY : scf.stats
scf.stats: $(SCF_STATS_TOCOMP)
endif

# individual rules
WRITE_STATS = scripts/write_stats.sh
# ABySS contigs or scaffolds;
%_ctg.stats %_scf.stats : %.fa
	$(WRITE_STATS) $< > $@

# SOAP contigs
%_ctg.stats : %.contig
	$(WRITE_STATS) $< > $@

# SOAP scaffolds
%_scf.stats : %.scafSeq
	$(WRITE_STATS) $< > $@
