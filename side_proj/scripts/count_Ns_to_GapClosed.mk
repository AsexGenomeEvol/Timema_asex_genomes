GENOMES := $(shell find data -name "*_GC.fa")
GENOME_Ns := $(patsubst %.fa, %_Ns.tsv, $(GENOMES))

.PHONY : count.all.GC.Ns
count.all.Ns : $(GENOME_Ns)

%_Ns.tsv : %.fa
	scripts/generic_genomics/fasta2number_of_N.sh $< > $@
