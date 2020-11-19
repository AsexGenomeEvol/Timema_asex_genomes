# author : Kamil S. Jaron
# ## are lines with help
# #! are lines that will be also printed for full help
## jobs for cluster uisng lsf scheduler are marked by symbol @ at the beginning of each help entry

include template.mk

## make help : print help
## make help.full : print complete help to all the automated steps including those that were not part of the final assembly or analyses

.PHONY : help.full
help.full : Makefile
	@sed -n 's/^#[#!]//p' $<

##
## A - organisation
## TODO make raw_reads : download raw reads from SRA
# data/%/checkpoints/raw_reads :
# 	download it > $@

##
## B - read parsing
## make trimmed_reads : @ trimm all reads (require raw_reads)
## make clean.raw_reads : delete raw reads that are already trimmed
## make stats.trimmed_reads : counts lengths of sequences in the trimmed read files
## make stats/reads/coverage.table.tsv : make a table of coverages of trimmed reads
## make kmergenie : @ predict otimal kmer using kmergenie from trimmed reads
#! make kmergenie.filt : @ predict otimal kmer using kmergenie from contamination-filtered trimmed reads
#! make stats/reads/kmertable_filt.tsv : make a table of optimal kmers from filtered reads
## make filtered_reads : make filtered reads (it filters reads that do not map to non-filtered scaffolds in version b3v6)
## make stats.filtered_reads : compute stats of filtered reads (a number of sequenced nts)
## make stats/reads/filtered_reads_coverage.table.tsv : make a table out of it (incl. proportion of surviving reads)
#! make rescue.reads : create filtered reads from bam files

.PHONY : trimmed_reads
.PHONY : clean.raw_reads
.PHONY : stats.trimmed_reads
.PHONY : stats/reads/coverage.table.tsv
.PHONY : kmergenie
.PHONY : kmergenie.filt
.PHONY : stats/reads/kmertable_filt.tsv

trimmed_reads clean.raw_reads stats.trimmed_reads stats/reads/coverage.table.tsv :
	$(MAKE) $@ -f B_read_parsing/trimming.mk

kmergenie kmergenie.filt stats/reads/kmertable_filt.tsv:
	$(MAKE) $@ -f B_read_parsing/kmergenie.mk

#stats/kmertable.tsv: $(KMERFILES)
#	touch stats/kmertable.tsv


##
## C - assembly
#! make assemble.batch1 : @ submit a batch of assemblies (SOAP, kmers 41, 43, 45, --nomse --nompe --fewdata)
#! make assemble.batch1.next_iter : @ do next iteration in otpimal asm search (SOAP, --nomse --nompe --fewdata)
#! make assemble.batch2.next_iter : @ do an iteration in otpimal asm search based on the optimal of b1 (SOAP, --nomse --nompe)
## make assemble.batch3 : @ submit a batch of assemblies (abyss, opt kmers, --nomse --nompe --fewdata)
## make clean.asm 2> stats/assemblies/current_optimal.tsv : clean all intermediate files from assmeblies that are suboptimal
#! make reassemble.batch3 : @ submit a batch of assemblies (abyss, opt kmers, --nomse --nompe --fewdata --filtered)

.PHONY : assemble.batch1
.PHONY : assemble.batch1.next_iter
.PHONY : assemble.batch2.next_iter
.PHONY : assemble.batch3
.PHONY : clean.asm
.PHONY : reassemble.batch3

assemble.batch1 assemble.batch1.next_iter \
assemble.batch2.next_iter \
assemble.batch3 reassemble.batch3 \
clean.asm :
	$(MAKE) $@ -f C_contig_assembly/assembly.mk

##
## D - scaffolding
#! make gapfilling.batch1 : @ gapfill batch 1 (2 is out of game)
## make filter.batch3 : @ filter contigs of batch 3 to greater than 250
## make indexed.batch3 : @ index filtered contigs of batch 3
## make mapping.batch3 : @ map pe reads to contig assembly of batch3
## make scaffold.batch3 : @ scaffold abyss contigs with BESST
## make gapfilling.batch3 : @ gapfill batch 3
#! make indexed.reasm.batch3 : @ index filtered contigs of reassembled batch 3
#! make mapping.reasm.batch3 : @ map pe reads to contig assembly of batch3
#! make scaffold.reasm.batch3 : @ scaffold reassembled abyss contigs with BESST
#! make gapfilling.reasm.batch3 : @ gapfill reassembled scaffolds

.PHONY : filter.batch3
.PHONY : indexed.batch3
.PHONY : indexed.reasm.batch3
.PHONY : mapping.batch3
.PHONY : mapping.reasm.batch3
.PHONY : scaffold.batch3
.PHONY : scaffold.reasm.batch3
.PHONY : gapfilling.batch1
.PHONY : gapfilling.batch3
.PHONY : gapfilling.reasm.batch3

filter.batch3 \
indexed.batch3 indexed.reasm.batch3 \
mapping.batch3 mapping.reasm.batch3 \
scaffold.batch3 scaffold.reasm.batch3 \
gapfilling.batch1 gapfilling.batch3 gapfilling.reasm.batch3 :
	$(MAKE) $@ -f D_scaffolding/scaffolding.mk

##
## E - evaluaiton
## make asm.stats  : will compute stats for all assemblies in data/*/assemblies
## make asm.tables : make all assembly tables (per species per ctg/scf)
#! make ctg.tables : make tables of assemblies of contigs per species
#! make scf.tables : make tables of assemblies of scaffolds per species
## make stats/assemblies/ctgs_fulltable.tsv -f assembly_benchmarking_stats.mk : make one big table of contig assemblies
## make stats/assemblies/scfs_fulltable.tsv -f assembly_benchmarking_stats.mk : make one big table of scaffold assemblies
# TODO nice benchmarking table (it got to be somewhere)
#! make filter.batch1 : @ filter scaffolds smaller than 250 bases from b1
## make busco : @ compute busco score for batch 1 and 3
## make quast : @ compute quast for batch 1 and 3
#! reapr STEP WAS NOT USED IN THE END (it had some troubles)
#! make reapr_map : @ map 5000 is library for the 1st and the 3rd batch
#! make reapr : @ call reapr on Gapfilled batch 1 and 3
## make counts.Ns : create tsv file with number of unknown nucleotides for each GapClosed genome
#! make stats/assemblies/batch_table.tsv
#! make archive.b1 : move b1 assemblies to archive (do not use this one)

.PHONY : asm.stats asm.tables ctg.tables scf.tables stats/assemblies/ctgs_fulltable.tsv stats/assemblies/scfs_fulltable.tsv
asm.stats asm.tables ctg.tables scf.tables stats/assemblies/ctgs_fulltable.tsv stats/assemblies/scfs_fulltable.tsv :
	$(MAKE) $@ -f E_assembly_evaluation/assembly_benchmarking_stats.mk

.PHONY : filter.batch1 filter.batch2 archive.b1
filter.batch1 filter.batch2 archive.b1 :
	$(MAKE) $@ -f E_assembly_evaluation/filtering.mk

.PHONY : busco
.PHONY : quast
busco quast :
	$(MAKE) $@ -s BATCH=1 -f E_assembly_evaluation/$@.mk
	$(MAKE) $@ -s BATCH=3 -f E_assembly_evaluation/$@.mk

.PHONY : reapr_map
.PHONY : reapr
reapr_map reapr :
	$(MAKE) reapr -s BATCH=1 SCRIPT=$@ -f E_assembly_evaluation/reapr.mk
	$(MAKE) reapr -s BATCH=3 SCRIPT=$@ -f E_assembly_evaluation/reapr.mk

.PHONY : counts.Ns
counts.Ns :
	$(MAKE) -f scripts/count_Ns_to_GapClosed.mk

stats/assemblies/batch_table.tsv : E_assembly_evaluation/make_batch_table.R
	Rscript $< $@


##
## F - reference polishing

## make rename.batch3 : rename headers of batch3 gapfilled assemblies
## make.reference.tab.b?v?? : creates a table of stats of a reference (stats/reference/*.tsv)
## make.reference.lengths.b?v?? : creates a table of stats of a reference (stats/reference/*.tsv)

.PHONY : rename.batch3
rename.batch3 :
	$(MAKE) $@ -f F_polishing_reference/renaming.mk

make.reference.tab.% :
	$(MAKE) stats/reference/$*.tsv -s VERSION=$* -f F_polishing_reference/reference_table.mk

make.reference.lengths.% :
	$(MAKE) lengths -s VERSION=$* -f F_polishing_reference/scaffold_lengths.mk

##
## G - mapping
## make map.ref.b?v?? : @ map pe libraries of reference individual to the specified version of asm
## make map.masked_ref.div?0 :  @ map pe libraries of reference individual to the specified masked version of b3v06 asm
## make map.to.contaminants : @ map pe libraries to contaminants found in b3v05
## make all.bam.stats : find all bam files in mapping subfolder and compute stats

map.ref.% :
	$(MAKE) map.all -s VERSION=$* IS=350 -f G_mapping/map_ref_reads_to_ref.mk
	$(MAKE) map.all -s VERSION=$* IS=550 -f G_mapping/map_ref_reads_to_ref.mk
	$(MAKE) map.all -s VERSION=$* IS=700 -f G_mapping/map_ref_reads_to_ref.mk

map.masked_ref.% :
	$(MAKE) map.all -s MASKED=$* IS=350 -f G_mapping/map_ref_reads_to_ref.mk
	$(MAKE) map.all -s MASKED=$* IS=550 -f G_mapping/map_ref_reads_to_ref.mk
	$(MAKE) map.all -s MASKED=$* IS=700 -f G_mapping/map_ref_reads_to_ref.mk

.PHONY : map.to.contaminants
map.to.contaminants :
	$(MAKE) map.to.contamination -s IS=350 -f G_mapping/map_ref_reads_to_ref.mk
	$(MAKE) map.to.contamination -s IS=550 -f G_mapping/map_ref_reads_to_ref.mk
	$(MAKE) map.to.contamination -s IS=700 -f G_mapping/map_ref_reads_to_ref.mk
	$(MAKE) map.mp.to.contamination -s IS=3000 -f G_mapping/map_ref_reads_to_ref.mk
	$(MAKE) map.mp.to.contamination -s IS=5000 -f G_mapping/map_ref_reads_to_ref.mk

.PHONY : all.bam.stats
all.bam.stats :
	$(MAKE) $@ -f G_mapping/stats_of_bam.mk


#clean.trimmed_reads.% : data/%/checkpoints/frozen data/%/checkpoints/mappings
#	rm all trimmed reads

##
## H - theta estimates
## make recalibrate : recalibrate quality scores of mapped div30 library
## make calc.thetas.div?0 : will calculate thetas for all bam files of filtred reads mapped to masked genome
## make calc.all_lib.recab.thetas  : will calculate thetas for merged bam files of filtred reads mapped to masked genome

.PHONY : recalibrate
.PHONY : calc.thetas.%
.PHONY : calc.all_lib.thetas.%

recalibrate :
	$(MAKE) call.all -f H_theta_calling/recalibrate.mk

calc.thetas.% :
	$(MAKE) call.all -s MASKED=$* -f H_theta_calling/call_atlas.mk

calc.all_lib.recab.thetas :
	$(MAKE) call.all -s ALL=yes -f H_theta_calling/call_atlas.mk


##
## I - theta analysis
## make plot.thetas.masked.div?0 :
## make plot.coverage.plots.div?0 :
## make heterozygosity.table :
## make plot.all.sp.coverage.plots.div?0 :
## make plot.violin.IND : where IND is an individual

.PHONY : plot.thetas.masked.%
plot.thetas.masked.% :
	$(MAKE) plot.thetas.masked -s MASKED=$* -f I_theta_analysis/analyze_thetas.mk

.PHONY : plot.coverage.plots.div?0
plot.coverage.plots.% :
	$(MAKE) plot.coverage.plots -s MASKED=$* -f I_theta_analysis/analyze_thetas.mk

.PHONY : heterozygosity.table
heterozygosity.table :
	$(MAKE) stats/b3v06_meidan_heterozygosity.tsv -f I_theta_analysis/analyze_thetas.mk

.PHONY : plot.all.sp.coverage.plots
plot.all.sp.coverage.plots.% :
	$(MAKE) plot.all.sp.coverage.plots -s MASKED=$* -f I_theta_analysis/analyze_thetas.mk

.PHONY : plot.violin.%
plot.violin.% :
	$(MAKE) plot.violin -s IND=$* -f I_theta_analysis/analyze_thetas.mk

##
## J - variant calling
## make call.gatk : call variants using GATK

.PHONY : call.gatk
call.gatk :
	$(MAKE) call.gatk -f J_SNP_calling/call.variants.mk

##
## N - genome alignment
## make genome.aligments : calculate genome alignments

.PHONY : genome.aligments
genome.aligments :
	$(MAKE) align_all -f N_genome_alignment/genome_alignment.mk

##
## scripts - indexing
## make index.all.bam : @ find all bam files
## make index.bwa.ref.b?v?? : @ will index a batch / versoion of all species (.bwa)
## make index.bwa.masked.div?0 : @ will index a batch of b3v06 softmasked genomes (.bwa)
## make clean.bwa.index.b?v?? : will clean bwa index of specified version (.bwa)
## make index.sam.ref.b?v?? : @ will index a batch / versoion of all species (.fai)

.PHONY : index.all.bam
index.all.bam :
	$(MAKE) -f scripts/index_all_bam.mk

index.bwa.ref.% :
	$(MAKE) index.all -s VERSION=$* -f scripts/index_bwa_reference.mk

index.bwa.masked.% :
	$(MAKE) index.all -s MASKED=$* -f scripts/index_bwa_reference.mk

clean.bwa.index.% :
	$(MAKE) clean.index -s VERSION=$* -f scripts/index_bwa_reference.mk

index.sam.ref.% :
	$(MAKE) index.all -s VERSION=$* -f scripts/index_sam_reference.mk
