# -*- coding: utf-8 -*-


#######################################################################
### Patrick Tran van : patrick.tranvan@unil.ch
###
### Pipeline for the assembly. See the github for usage.
###
#######################################################################

import os
import re
import argparse
import sys
import glob
import os.path
import random
from Bio import SeqIO

#### HOW TO USE

## python $SC/asv3.py -s fastuniq -i1 raw_read_merged/ -o fastuniq.sh
## python $SC/asv3.py -s nxtrim -i1 raw_reads/ -e patrick.tranvan@unil.ch -o nxtrim.sh
## python $SC/asv3.py -s cat_nxtrim -i1 nxtrim/ -e patrick.tranvan@unil.ch -o cat_nxtrim.sh
## python $SC/asv3.py -s trimmomatic -i1 raw_reads/ -i2 merge/ -e patrick.tranvan@unil.ch -o trimmomatic.sh
## python $SC/asv3.py -s MP_stat -i1 nxtrim/ -o MP.stat.txt
## python $SC/asv3.py -s PE_stat -i1 trim/ -o PE.stat.txt
## python $SC/asv3.py -s taxonomy_list -i1 Sm2.blobDB.table.txt
## python $SC/asv3.py -s blobtool_contaminants -i1 Sm2.blobDB.table.txt
## python $SC/asv3.py -s link_contaminants -i1 specie_contaminant_select_unique.txt -i2 fungi_bacteria_select.txt
## python $SC/asv3.py -s contaminants_filtration -i1 contaminants.fasta -o contaminants_fitered.fasta
## python $SC/asv3.py -s bwa -i1 reads.list -i2 db -e my_email@unil.ch -o bwa.sh
## python $SC/asv3.py -s bbmap -i1 reads.list -e patrick.tranvan@unil.ch -o bbmap.sh
## python $SC/asv3.py -s bbmap_sam -i1 reads.list -e patrick.tranvan@unil.ch -o bbmap_sam.sh
## python $SC/assembly.py -s samtofastq -i1 bowtie2/ -e my_email@unil.ch -o samtofastq.sh
## python $SC/asv3.py -s sort -i1 bwa/ -e patrick.tranvan@unil.ch -o sort.sh
## python $SC/asv3.py -s sort_query -i1 bwa/ -e patrick.tranvan@unil.ch -o sort_query.sh
## python $SC/asv3.py -s flagstat -i1 bwa/ -e patrick.tranvan@unil.ch -o flagstat.sh
## python $SC/asv3.py -s map_filter -i1 bwa/ -e patrick.tranvan@unil.ch -o map_filter.sh
## python $SC/asv3.py -s bamtofastq -i1 bwa/ -e patrick.tranvan@unil.ch -o bamtofastq.sh
## python $SC/asv3.py -s fastq_stat -i1 reads/ -e patrick.tranvan@unil.ch -o fastq_stat.sh
## python $SC/asv3.py -s fastq_stat_table -i1 reads/ -o fastq.stat.txt
## python $SC/asv3.py -s flag_fastq_stat -i1 <bwa_directory> -i2 <read_directory> -o flag_fastq_stat.txt
## python $SC/asv3.py -s mark_duplicates -i1 <bwa_directory> -e patrick.tranvan@unil.ch -o mark_duplicates.sh
## python $SC/asv3.py -s rename -i1 fasta -o fasta_rename.fasta
## python $SC/asv3.py -s rename_TE -i1 fasta -o fasta_rename.fasta
## python $SC/asv3.py -s reformat -i1 fasta -o fasta_rename.fasta
## python $SC/asv3.py -s busco_assignment -i1 full_table_output.tsv -i2 tax_list.txt -o busco_assignment.txt
## python $SC/asv3.py -s cegma_assignment -i1 1_Tps.cegma.gff -i2 tax_list.txt -o cegma_assignment.txt
## python $SC/asv3.py -s busco_coverage -i1 busco_tsv -i2 coverage_file -i3 fasta_file -o busco_coverage.txt
## python $SC/asv3.py -s cegma_reformat -i1 output.cegma.id -i2 output.cegma.gff -o cegma_summary.txt
## python $SC/asv3.py -s cegma_coverage -i1 cegma_summary -i2 coverage_file -i3 fasta_file -o cegma_coverage.txt
## python $SC/asv3.py -s diff_scaf_busco_assign -i1 scaf_file -i2 busco_file 
## python $SC/asv3.py -s busco_cegma_cov -i1 busco_coverage.txt -i2 cegma_coverage.txt -o busco_cegma_cov.txt
## python $SC/asv3.py -s random_coverage -i1 fasta_file -o output_file
## python $SC/asv3.py -s scaffolds_busco_cegma -i1 busco_tsv -i2 cegma.gff -i3 contaminant_scaffolds -o scaffold_list.txt
## python $SC/asv3.py -s cat_reseq -i1 reads -e patrick.tranvan@unil.ch -o cat.sh
## python $SC/asv3.py -s bwa_reseq -i1 trim/ -i2 db -i3 Sm -o bwa_reseq.sh
## python $SC/asv3.py -s coverage_manual -i1 stat_pasrse -i2 contaminant_scaffolds.txt -o stat_parse.txt
## python $SC/asv3.py -s blast_manual -i1 blast -i2 contaminant_scaffolds.txt -o .vs.nt.cul5.1e25.megablast.out
## python $SC/asv3.py -s read_size -i1 pwd -o read_size.sh

def fastuniq(raw_read_merged_directory, output_file):

	"""
	Write a script to run FastUniq.

	1) The path has to be like this:
	raw_reads/<ins_length>/*.fastq*
	
	2) Name of PE reads must be the same except the end:
	repo/3000/As6.raw.550.pair1.fastq 
	repo/3000/As6.raw.550.pair2.fastq
		
	"""
		
	if os.path.isfile(output_file):
		os.remove(output_file)
	
	os.makedirs("fastuniq/180/")
	os.makedirs("fastuniq/350/")
	os.makedirs("fastuniq/550/")
	os.makedirs("fastuniq/3000/")
		
	script = open(output_file,"a")
	script.write("#!/bin/bash\n\n")
	script.write("source /scratch/beegfs/monthly/ptranvan/Software/FastUniq/1.1.sh\n\n")
	script.close()
	

	PE_data_repo = os.path.join(raw_read_merged_directory, "*")
	
	i = 0
	while i<len(glob.glob(PE_data_repo)):
	
		ins_length_path = sorted(glob.glob(PE_data_repo))[i]	
		ins_length = ins_length_path.split("/")[-1]
					
		read_repo = os.path.join(ins_length_path, "*pair1*")
		

		j = 0
		while j<len(glob.glob(read_repo)):
		
		
			pe_data = sorted(glob.glob(read_repo))[j]
			sample_name = (pe_data.split("/")[-1]).split(".")[0]
			sample_name_insert = (pe_data.split("/")[-1]).split(".pair1.fastq")[0]
			#print sample_name	
			fastuniq_ins_length_path = os.path.join("fastuniq", ins_length)		
			#print pe_data


			script = open(os.path.join(fastuniq_ins_length_path, sample_name_insert + ".input"), "a")			
			script.write("{}\n{}".format(pe_data, pe_data.replace("pair1", "pair2")))					
			script.close()

		
			script = open(output_file,"a")
			
			script.write("echo '{}'\nfastuniq -i {} -t q -o {} -p {} -c 0\n\n".format(sample_name_insert, os.path.join(fastuniq_ins_length_path, sample_name_insert + ".input"), os.path.join(fastuniq_ins_length_path, "{}.raw.nodup.{}.pair1.fastq".format(sample_name, ins_length)), os.path.join(fastuniq_ins_length_path, "{}.raw.nodup.{}.pair2.fastq".format(sample_name, ins_length))))
					
			script.close()

			j+=1

		i += 1
		
				
	print "Script has been written."
				 
def nxtrim(raw_read_directory, email, output_file):

	"""	
	Write a script to run nextclip in Vital-IT cluster.
	
	1) The path has to be like this:
	raw_reads/<ins_length>/*.fastq*
	
	2) Name of PE reads must be the same except the end:
	repo/3000/As6.raw.3000.pair1.fastq 
	repo/3000/As6.raw.3000.pair2.fastq
	"""
		
	if os.path.isfile(output_file):
		os.remove(output_file)
		
	script = open(output_file,"a")
	script.write("#!/bin/bash\n\n")
	script.close()
	
	data_repo = os.path.join(raw_read_directory, "*")
	
	i = 0
	while i<len(glob.glob(data_repo)):
	
		ins_length_path = sorted(glob.glob(data_repo))[i]	
		ins_length = ins_length_path.split("/")[-1]
		
		if int(ins_length) > 1000:	# Works only for mate pair
		
			read_repo = os.path.join(ins_length_path, "*.pair1.fastq")
			
			j = 0
			while j<len(glob.glob(read_repo)):
			
				mp_data = sorted(glob.glob(read_repo))[j]
				sample_name = (mp_data.split("/")[-1]).split(".")[0]
				print sample_name	
				nxtrim_ins_length_path = os.path.join("nxtrim", ins_length)		
				#print mp_data
				
				
				script = open(output_file,"a")
								
				#script.write("echo 'mkdir -p {} && module add UHTS/Quality_control/NxTrim/0.4.1 && nxtrim -1 {} -2 {} -O {} --separate' | bsub -q dee-hugemem -M 10000000 -R 'span[ptile=10] rusage[mem=10000]' -u {} -N -J nxtrim_{} -e {}\n\n".format(nxtrim_ins_length_path, mp_data, mp_data.replace("pair1", "pair2"), os.path.join(nxtrim_ins_length_path, sample_name), email, sample_name, os.path.join(nxtrim_ins_length_path, sample_name + ".nxtrim.log")))
				
				script.write("echo 'mkdir -p {} && module add UHTS/Quality_control/NxTrim/0.4.1 && nxtrim -1 {} -2 {} -O {} --separate --preserve-mp --minlength 40' | bsub -q dee-hugemem -M 10000000 -R 'span[ptile=10] rusage[mem=10000]' -u {} -N -J nxtrim_{} -e {}\n\n".format(nxtrim_ins_length_path, mp_data, mp_data.replace("pair1", "pair2"), os.path.join(nxtrim_ins_length_path, sample_name), email, sample_name, os.path.join(nxtrim_ins_length_path, sample_name + ".nxtrim.log")))
												
				script.close()
				
				j+=1
			
		i += 1
				
	print "Script has been written."

def cat_nxtrim(nxtrim_directory, email, output_file):

	"""	
	Concatenate categories mp and unknow (real mate pair) (FR oriented) after nxtrim jobs.
	"""
		
	if os.path.isfile(output_file):
		os.remove(output_file)
		
	script = open(output_file,"a")
	script.write("#!/bin/bash\n\n")
	script.close()
	
	data_repo = os.path.join(nxtrim_directory, "*")
	
	i = 0
	while i<len(glob.glob(data_repo)):
	
		ins_length_path = sorted(glob.glob(data_repo))[i]	
		ins_length = ins_length_path.split("/")[-1]
				
		read_repo = os.path.join(ins_length_path, "*R1.mp.fastq.gz")
		
		j = 0
		while j<len(glob.glob(read_repo)):
		
			mp_data = sorted(glob.glob(read_repo))[j]
			sample_name = (mp_data.split("/")[-1]).split(".")[0][:-3]
			#print sample_name	
			nxtrim_ins_length_path = os.path.join(nxtrim_directory, ins_length)		
			#print mp_data
			
			script = open(output_file,"a")
			
			script.write("echo 'mkdir -p {} && zcat {} {} > {}' | bsub -u {} -N -J cat_{}\n\n".format(os.path.join(nxtrim_ins_length_path, "merge"), mp_data, mp_data.replace("R1.mp.fastq.gz", "R1.unknown.fastq.gz"), os.path.join(nxtrim_ins_length_path, os.path.join("merge", sample_name + ".nxtrim.3000.pair1.fastq")), email, sample_name))

			script.write("echo 'mkdir -p {} && zcat {} {} > {}' | bsub -u {} -N -J cat_{}\n\n".format(os.path.join(nxtrim_ins_length_path, "merge"), mp_data.replace("R1.mp.fastq.gz", "R2.mp.fastq.gz"), mp_data.replace("R1.mp.fastq.gz", "R2.unknown.fastq.gz"), os.path.join(nxtrim_ins_length_path, os.path.join("merge", sample_name + ".nxtrim.3000.pair2.fastq")), email, sample_name))

			script.close()

			j += 1
			
		i += 1
				
	print "Script has been written."
		
def trimmomatic(raw_read_directory, nxtrim_merge_directory, email, output_file):

	"""
	Write a script to run trimmomatic in Vital-IT cluster.

	1) The path has to be like this:
	raw_reads/<ins_length>/*.fastq*
	
	2) Name of PE reads must be the same except the end:
	repo/3000/As6.raw.550.pair1.fastq 
	repo/3000/As6.raw.550.pair2.fastq
		
	"""
		
	if os.path.isfile(output_file):
		os.remove(output_file)
		
	script = open(output_file,"a")
	script.write("#!/bin/bash\n\n")
	script.close()
	
	# Job for PE
	PE_data_repo = os.path.join(raw_read_directory, "*")
	
	i = 0
	while i<len(glob.glob(PE_data_repo)):
	
		ins_length_path = sorted(glob.glob(PE_data_repo))[i]	
		ins_length = ins_length_path.split("/")[-1]
		
		if int(ins_length) < 1000:
			
			read_repo = os.path.join(ins_length_path, "*pair1*")
			

			j = 0
			while j<len(glob.glob(read_repo)):
			
			
				pe_data = sorted(glob.glob(read_repo))[j]
				sample_name = (pe_data.split("/")[-1]).split(".")[0]
				#print sample_name	
				trim_ins_length_path = os.path.join("trim", ins_length)		
				#print pe_data
				
				# rnaseq parameter
				
				script = open(output_file,"a")
				
				script.write("echo 'mkdir -p {} && module add Development/java/1.8.0_152 && java -jar /scratch/beegfs/monthly/ptranvan/Software/trimmomatic/0.36/trimmomatic.jar PE -threads 1 -phred33 {} {} {} {} {} {} ILLUMINACLIP:/scratch/beegfs/monthly/ptranvan/Software/trimmomatic/0.36/adapters/all-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:3:20 MINLEN:80' | bsub -q dee-hugemem -M 20000000 -R 'span[ptile=1] rusage[mem=20000]' -u {} -N -J trimmomatic_{} -e {}\n\n".format(trim_ins_length_path, pe_data, pe_data.replace("pair1", "pair2"), os.path.join(trim_ins_length_path, "{}.trim.pair1.fastq.gz".format(sample_name)), os.path.join(trim_ins_length_path, "{}.trim.single1.fastq.gz".format(sample_name)), os.path.join(trim_ins_length_path, "{}.trim.pair2.fastq.gz".format(sample_name)), os.path.join(trim_ins_length_path, "{}.trim.single2.fastq.gz".format(sample_name)), email, sample_name, os.path.join(trim_ins_length_path, sample_name + ".trimmomatic.log")))
						
				script.close()

				j+=1

		i += 1
	
	# Job for MP
	"""
	MP_data_repo = os.path.join(nxtrim_merge_directory, "*pair1.fastq")
	
	i = 0
	while i<len(glob.glob(MP_data_repo)):
			
			
		mp_data = sorted(glob.glob(MP_data_repo))[i]
		sample_name = (mp_data.split("/")[-1]).split(".")[0]
		#print sample_name	

		
		script = open(output_file,"a")
		
		script.write("echo 'java -jar /scratch/beegfs/monthly/ptranvan/Software/trimmomatic/0.36/trimmomatic.jar PE -threads 1 -phred33 {} {} {} {} {} {} ILLUMINACLIP:/scratch/beegfs/monthly/ptranvan/Software/trimmomatic/0.36/adapters/all-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:60' | bsub -q dee-hugemem -M 30000000 -R 'span[ptile=1] rusage[mem=30000]' -u {} -N -J trimmomatic_{} -e {}\n\n".format(mp_data, mp_data.replace("pair1.fastq", "pair2.fastq"), sample_name + ".nxtrim.trim.3000.pair1.fastq", sample_name + ".nxtrim.trim.3000.single1.fastq", sample_name + ".nxtrim.trim.3000.pair2.fastq", sample_name + ".nxtrim.trim.3000.single2.fastq", email, sample_name, sample_name + ".trimmomatic.log"))
				
		script.close()
		
		i+=1	
	
	"""			
	print "Script has been written."	
	
def MP_stat(nxtrim_directory, output_file):

	"""
	Statistics for Mate-Pair, after nxtrim()

	"""
	
	if os.path.isfile(output_file):
		os.remove(output_file)
		
	data_repo = os.path.join(nxtrim_directory, "*")
	
	script = open(output_file,"a")
	script.write("Lib\tSample\tInput read pairs\tMP\tPE\tUnknown\tSingle\n")
	
	i = 0
	while i<len(glob.glob(data_repo)):
	
		ins_length_path = sorted(glob.glob(data_repo))[i]	
		read_repo = os.path.join(ins_length_path, "*nxtrim.log")
		
		j = 0
		while j<len(glob.glob(read_repo)):
		
			stat = []
			nxtrim_data = sorted(glob.glob(read_repo))[j]
			
			insert_size = nxtrim_data.split("/")[-2]
			#print insert_size
			sample_name = (nxtrim_data.split("/")[-1]).split(".")[0]
			#print sample_name
			print nxtrim_data
			
			nxtrim_output_log = open(nxtrim_data,"r")
			#stat = ([line.split() for line in nxtrim_output_log if re.search("Trimming summary",line)])
			#print stat
			for line in nxtrim_output_log:
				if  re.search("Trimming summary",line):
				
					stat.append(next(nxtrim_output_log).split()[0])
					next(nxtrim_output_log)
					next(nxtrim_output_log)
					next(nxtrim_output_log)
					next(nxtrim_output_log)
					
					current_stat = next(nxtrim_output_log)
					
					stat.append(current_stat.split()[0])
					stat.append(current_stat.split()[4])
					
					current_stat = next(nxtrim_output_log)
					
					stat.append(current_stat.split()[0])
					stat.append(current_stat.split()[4])
					
					current_stat = next(nxtrim_output_log)
					
					stat.append(current_stat.split()[0])
					stat.append(current_stat.split()[4])
					
					current_stat = next(nxtrim_output_log)
					
					stat.append(current_stat.split()[0])
					stat.append(current_stat.split()[4])
					
			nxtrim_output_log.close()
			
			#print stat
			
			
			script.write("{}\t{}\t{}\t{} ({})\t{} ({})\t{} ({})\t{} ({})\n".format(insert_size, sample_name, int(stat[0]), int(stat[1]), stat[2], int(stat[3]), stat[4], int(stat[5]), stat[6], int(stat[7]), stat[8]))
			
			j+=1
			
		i += 1
	
	script.close()
	
	print "Table has been written."

def PE_stat(trim_directory, output_file):

	"""
	Statistics for Paired-end after trimmomatic().
	"""
	
	if os.path.isfile(output_file):
		os.remove(output_file)
		
	data_repo = os.path.join(trim_directory, "*")
	
	script = open(output_file,"a")
	script.write("Lib\tSample\tInput read pairs\tBoth surviving\tF only surviving\tR only surviving\tDropped\n")
	
	i = 0
	while i<len(glob.glob(data_repo)):
	
		ins_length_path = sorted(glob.glob(data_repo))[i]	
		read_repo = os.path.join(ins_length_path, "*trimmomatic.log")
		
		j = 0
		while j<len(glob.glob(read_repo)):
		
			trim_data = sorted(glob.glob(read_repo))[j]
			#print pe_data # trim/350/160603_SND393_A_L002_HYI-120_output_log.txt
			
			insert_size = trim_data.split("/")[-2]
			#print insert_size
			
			
			sample_name = (trim_data.split("/")[-1]).split(".")[0]
			#print sample_name
			
			trimmomatic_output_log = open(trim_data,"r")
			stat = ([line.split() for line in trimmomatic_output_log if re.search("Input Read Pairs",line)][0])
			trimmomatic_output_log.close()
			
			script.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(insert_size, sample_name, stat[3], stat[6] + " " + stat[7], stat[11] + " " + stat[12], stat[16] + " " + stat[17], stat[19] + " " + stat[20]))
		
			j+=1
			
		i += 1
	
	script.close()

	print "Table has been written."
	
def taxonomy_list(blob_table):
	"""
	Create a table with the scaffold name and it's taxonomic assignation.
	"""

	if os.path.isfile("tax_list.txt"):
		os.remove("tax_list.txt")
			
	not_contaminants_domain = ["Eukaryota", "no-hit"]
            		
	blob_table_file = open(blob_table,"r")
	for scf in blob_table_file:
	
		if '#' in scf:
			pass
		else:
			scaffold = re.split("[\r\t\n]",scf)[0]
			domain = re.split("[\r\t\n]",scf)[5]
			phylum = re.split("[\r\t\n]",scf)[9]
			specie = re.split("[\r\t\n]",scf)[25]
			
			
			#print scf.split()
	
			script3 = open("tax_list.txt","a")			
			script3.write("{}\t{}\n".format(scaffold, phylum))
			script3.close()
			
	blob_table_file.close()

	print "\n\nTable has been written.\n"
			
def blobtool_contaminants(blob_table):
	"""
	Make a list of contaminants scaffolds for seqtk subseq.
	"""
	
	specie_contaminant_wide = []
	specie_contaminant_wide_dic = {}
            		
	blob_table_file = open(blob_table,"r")
	for scf in blob_table_file:
	
		if '#' in scf:
			pass
		else:
			scaffold = re.split("[\r\t\n]",scf)[0]
			domain = re.split("[\r\t\n]",scf)[5]
			phylum = re.split("[\r\t\n]",scf)[9]
			phylum_hits = re.split("[\r\t\n]",scf)[12]
			# tax0=Streptophyta:837.0|Chordata:561.0|Eukaryota-undef:182.0|Ascomycota:165.0;

			order = re.split("[\r\t\n]",scf)[13]
			genus = re.split("[\r\t\n]",scf)[21]
			specie = re.split("[\r\t\n]",scf)[25]
			
			#sinergia + thrips
			
			not_contaminants_phylums_wide = ["Annelida", "Arthropoda", "Brachiopoda", "Chordata", "Cnidaria", "Echinodermata", "Hemichordata", "Mollusca", "Nematoda", "Nemertea", "no-hit", "Onychophora", "Placozoa", "Platyhelminthes", "Porifera", "Priapulida", "Rotifera", "Tardigrada"]	#sm2
			
			# fsel
			
			#not_contaminants_phylums_wide = ["Arthropoda", "no-hit"]
				
			## Wide
			
			multi_hit_wide = False
			
			for not_contaminants_phylum_wide in not_contaminants_phylums_wide:
				if not_contaminants_phylum_wide in phylum_hits:
					multi_hit_wide = True	#hit include not_contaminants_phylums_wide
					break
					
			if multi_hit_wide == False:
			
				script = open("contaminant_scaffolds.txt","a")			
				script.write("{}\n".format(scaffold))
				script.close()
				
				script = open("contaminant_tax_list.txt","a")			
				script.write("{}\t{}\t{}\t{}\n".format(scaffold, phylum, genus, specie))
				script.close()
				
				specie_contaminant_wide.append(specie)
				specie_contaminant_wide_dic[specie] = phylum
			
	blob_table_file.close()

	for specie in list(set(specie_contaminant_wide)):
	
		script = open("unique_contaminant_specie_list.txt","a")			
		script.write("{}\n".format(specie))
		script.close()
		
		script = open("unique_contaminant_specie_phylum_list.txt","a")			
		script.write("{}\t{}\n".format(specie, specie_contaminant_wide_dic[specie]))
		script.close()
	
	print "\n\nTable has been written.\n"


def link_contaminants(specie_contaminant_unique_file, fungi_bacteria_select_file):
	"""
	Link to download contaminants genome.
	If complete genome available, download only complete genome else download everything.
	"""
	
	specie_contaminant_uniques = open(specie_contaminant_unique_file,"r")
	
	for specie_contaminant_unique in specie_contaminant_uniques:
	
		specie_contaminant = re.split("[\r\t\n]",specie_contaminant_unique)[0]	#name of the specie
		
		fungi_bacteria_selects = open(fungi_bacteria_select_file,"r")
		
		status_list = []
		link_list = []
		exist = False
		
		for fungi_bacteria_select in fungi_bacteria_selects:
		
			specie = re.split("[\r\t\n]",fungi_bacteria_select)[7]
			status = re.split("[\r\t\n]",fungi_bacteria_select)[11]	# complete genome or not
			link = re.split("[\r\t\n]",fungi_bacteria_select)[19]
			
			if specie_contaminant in specie:
				status_list.append(status)
				link_list.append(link)
				exist = True	# if the specie exists in the genbank database
			
		fungi_bacteria_selects.close()
		
		if exist:
		
			if status_list.count("Complete Genome") > 0:
				print specie_contaminant
				index_complete_genome = [i_complete for i_complete, status_data in enumerate(status_list) if status_data == "Complete Genome"]	# index list of complete genome
				print index_complete_genome
				
				for index in index_complete_genome:
					# Download all the complete genome that could correspond to different strains.
			
					script = open("link.txt","a")			
					script.write("{}/{}_genomic.fna.gz\n".format(link_list[index], link_list[index].split("/")[-1]))
					script.close()
					
			else:
			
				for index_no_complete in link_list:
				
					script = open("link.txt","a")			
					script.write("{}/{}_genomic.fna.gz\n".format(index_no_complete, index_no_complete.split("/")[-1]))
					script.close()		
	
	specie_contaminant_uniques.close()
	
	print "Done."

def contaminants_filtration(contaminants_file, contaminants_filtered_file):
	"""
	Remove reads with len = 0
	"""
	
	seq_record_sample=[]

	for seq_record in SeqIO.parse(contaminants_file, "fasta"):
		if len(seq_record.seq)== 0:
			print seq_record
		else:
			seq_record_sample.append(seq_record)
						
	SeqIO.write(seq_record_sample, contaminants_filtered_file, "fasta")
	

def bwa(reads_list, db, email, output_file):
	"""
	Map reads 
	reads_list: Reads path per line.
	No line at the end
	
	Different jobs between PE and single.
	"""
	
	output_dir = "bwa"
	
	if os.path.isfile(output_file):
		os.remove(output_file)
	
	script = open(output_file,"a")
	script.write("#!/bin/bash\n\n")
	script.close()
	
	single = ["all"]
	pair = ["pair1", "pair2"]
	
	trim_reads_file = open(reads_list,"r")
	content = trim_reads_file.readlines()

	i = 0 
	while i < len(content):
	
		sample_name = content[i].split("/")[-1].split(".fastq")[0]
		sample_type = content[i].split("/")[-1].split(".")[-2]
		#print sample_name
		#print sample_type
		#print content[i].split()
		
		if sample_type in single:
			
			script = open(output_file,"a")						
			script.write("echo 'mkdir -p {} && module add UHTS/Aligner/bwa/0.7.13 && bwa mem -M {} {} > {}' | bsub -q dee-hugemem -M 8000000 -R 'span[ptile=8] rusage[mem=8000]' -u {} -N -J bwa_{} -e {}\n\n".format(output_dir, db, content[i].split()[0], os.path.join(output_dir, sample_name + ".sam"), email, sample_name, os.path.join(output_dir, sample_name + ".bwa.output.log")))		
			script.close()
			
			i += 1
			
		else:	#PE
		
			#### Attention: a revoir la commande car normalement je veux du .sam
			script = open(output_file,"a")						
			script.write("echo 'mkdir -p {} && module add UHTS/Aligner/bwa/0.7.13 && module add UHTS/Analysis/samtools/1.3 && bwa mem -M {} {} {} -t 15 | samtools view -bS - > {}' | bsub -q dee-hugemem -M 10000000 -R 'span[ptile=15] rusage[mem=10000]' -n 15 -u {} -N -J bwa_{} -e {}\n\n".format(output_dir, db, content[i].split()[0], (content[i].split()[0]).replace("pair1", "pair2"), os.path.join(output_dir, sample_name + "2.bam"), email, sample_name, os.path.join(output_dir, sample_name + "2.bwa.output.log")))		
			script.close()			
			
			i += 2
			
	trim_reads_file.close()


	print "\nScript has been written.\n"

def bbmap(reads_list, email, output_file):
	"""
	Map reads with bbmap
	reads_list: Reads path per line.
	No line at the end
	
	Different jobs between PE and single.
	"""
	
	output_dir = "bwa"
	
	if os.path.isfile(output_file):
		os.remove(output_file)
	
	script = open(output_file,"a")
	script.write("#!/bin/bash\n\n")
	script.close()
	
	pair = ["pair1", "pair2"]
	
	trim_reads_file = open(reads_list,"r")
	content = trim_reads_file.readlines()

	i = 0 
	while i < len(content):
	
		file_name = content[i].split("/")[-1]
		sample_name = ".".join((file_name.split(".pair1")[0].split("."))[0:-1])
		insert_size = (file_name.split(".pair1")[0].split("."))[-1]
		#print sample_name
		#print sample_type
		#print content[i].split()

		
		script = open(output_file,"a")	
							
		script.write("echo 'module add UHTS/Aligner/BBMap/36.59 && bbmap.sh in={} in2={} outm1={} outm2={} outu={} killbadpairs=t pairedonly=t ihist={} -Xmx30g threads=30' | bsub -q dee-hugemem -M 30000000 -R 'span[ptile=32] rusage[mem=30000]' -n 32 -u {} -N -J bbmap_{} -e {}\n\n".format(content[i].split()[0], (content[i].split()[0]).replace("pair1", "pair2"), sample_name + ".filtered." + insert_size + ".pair1.fastq", sample_name + ".filtered." + insert_size + ".pair2.fastq", sample_name + ".contaminant." + insert_size + ".fastq", sample_name + insert_size + ".hist.txt", email, sample_name + insert_size, os.path.join(sample_name + insert_size + ".bbmap.log")))	
			
		script.close()			
		
		i += 2
			
	trim_reads_file.close()

	print "\nScript has been written.\n"

def bbmap_sam(reads_list, email, output_file):
	"""
	Map reads with bbmap to output sam files
	reads_list: Reads path per line.
	No line at the end
	
	"""
	
	output_dir = "bbmap"
	
	if os.path.isfile(output_file):
		os.remove(output_file)
	
	script = open(output_file,"a")
	script.write("#!/bin/bash\n\n")
	script.close()
		
	trim_reads_file = open(reads_list,"r")
	content = trim_reads_file.readlines()

	i = 0 
	while i < len(content):
	
		file_name = content[i].split("/")[-1]
		
		if file_name.split(".")[-2] == "pair12":
			
			sample_name = ".".join((file_name.split(".pair12")[0].split("."))[0:-1])
			insert_size = (file_name.split(".pair12")[0].split("."))[-1]
			#print sample_name
			#print sample_type
			#print content[i].split()
	
			
			script = open(output_file,"a")	
								
			script.write("echo 'module add UHTS/Aligner/BBMap/36.59 && bbmap.sh in={} outm={} outu={} requirecorrectstrand=f covstats={} -Xmx25g threads=30 delcoverage=f' | bsub -q dee-hugemem -M 25000000 -R 'span[ptile=32] rusage[mem=25000]' -n 30 -u {} -N -J bbmap_{} -e {}\n\n".format(content[i].split()[0], sample_name + ".mapped." + insert_size + ".sam", sample_name + ".unmapped." + insert_size + ".sam", sample_name + insert_size + ".cov.txt", email, sample_name + insert_size, os.path.join(sample_name + insert_size + ".bbmap_sam.log")))	
				
			script.close()			
			
			i += 1
			
		else:
		
			sample_name = ".".join((file_name.split(".pair1")[0].split("."))[0:-1])
			insert_size = (file_name.split(".pair1")[0].split("."))[-1]
			#print sample_name
			#print sample_type
			#print content[i].split()
	
			
			script = open(output_file,"a")	
								
			#script.write("echo 'module add UHTS/Aligner/BBMap/36.59 && bbmap.sh in={} in2={} outm={} outu={} requirecorrectstrand=f covstats={} -Xmx25g threads=30 delcoverage=f' | bsub -q dee-hugemem -M 25000000 -R 'span[ptile=32] rusage[mem=25000]' -n 32 -u {} -N -J bbmap_{} -e {}\n\n".format(content[i].split()[0], (content[i].split()[0]).replace("pair1", "pair2"), sample_name + ".mapped." + insert_size + ".sam", sample_name + ".unmapped." + insert_size + ".sam", sample_name + insert_size + ".cov.txt", email, sample_name + insert_size, os.path.join(sample_name + insert_size + ".bbmap_sam.log")))	
			
			##2

			script.write("echo 'mkdir -p {} && module add UHTS/Aligner/BBMap/36.59 && bbmap.sh in={} in2={} outm={} outu={} covstats={} killbadpairs=t pairedonly=t -Xmx25g threads=30 delcoverage=f' | bsub -q dee-hugemem -M 25000000 -R 'span[ptile=32] rusage[mem=25000]' -n 32 -u {} -N -J bbmap_{} -e {}\n\n".format(output_dir, content[i].split()[0], (content[i].split()[0]).replace("pair1", "pair2"), os.path.join(output_dir, sample_name + ".mapped." + insert_size + ".sam"), os.path.join(output_dir, sample_name + ".unmapped." + insert_size + ".sam"), os.path.join(output_dir, sample_name + insert_size + ".cov.txt"), email, sample_name + insert_size, os.path.join(output_dir, sample_name + insert_size + ".bbmap_sam.log")))	
							
			script.close()			
			
			i += 2
			
	trim_reads_file.close()

	print "\nScript has been written.\n"
	
		
def sort(mapping_directory, email, output_file):
	"""
	Using Samtools to sort bam file by leftmost coordinates (default).
	"""
	
	if os.path.isfile(output_file):
		os.remove(output_file)
	
	script = open(output_file,"a")
	script.write("#!/bin/bash\n\n")
	script.close()

	# For PE
	mapping_repo = os.path.join(mapping_directory, "*.bam")
	sort_path = os.path.join(mapping_directory, "sort")
	
	i = 0
	while i<len(glob.glob(mapping_repo)):

		sample = sorted(glob.glob(mapping_repo))[i]	
		file_name = sample.split("/")[-1]
		sample_name = ".".join((file_name.split(".pair12")[0].split("."))[0:-1])
		insert_size = (file_name.split(".pair12")[0].split("."))[-1]

		#print sample
			
		script = open(output_file,"a")
		
		script.write("echo 'mkdir -p {} && module add UHTS/Analysis/samtools/1.3 && samtools sort {} -o {}' | bsub -q dee-hugemem -M 20000000 -R 'span[ptile=1] rusage[mem=20000]' -u {} -N -J sort_{}\n\n".format(sort_path, sample, os.path.join(sort_path, sample_name + ".sorted." + insert_size + ".pair12.bam"), email, sample_name))		
	
		script.close()
		
		i += 1

	print "\nScript has been written.\n"

def sort_query(mapping_directory, email, output_file):
	"""
	Using Samtools to sort bam file by query name.
	"""
	
	if os.path.isfile(output_file):
		os.remove(output_file)
	
	script = open(output_file,"a")
	script.write("#!/bin/bash\n\n")
	script.close()

	# For PE
	mapping_repo = os.path.join(mapping_directory, "*.bam")
	sort_path = os.path.join(mapping_directory, "sort_query")
	
	i = 0
	while i<len(glob.glob(mapping_repo)):

		sample = sorted(glob.glob(mapping_repo))[i]	
		file_name = sample.split("/")[-1]
		sample_name = ".".join((file_name.split(".pair12")[0].split("."))[0:-1])
		insert_size = (file_name.split(".pair12")[0].split("."))[-1]

		#print sample
			
		script = open(output_file,"a")
		
		script.write("echo 'mkdir -p {} && module add UHTS/Analysis/samtools/1.3 && samtools sort -n {} -o {}' | bsub -q dee-hugemem -M 5000000 -R 'span[ptile=1] rusage[mem=5000]' -u {} -N -J sort_query_{}\n\n".format(sort_path, sample, os.path.join(sort_path, sample_name + ".sorted_query." + insert_size + ".pair12.bam"), email, sample_name))		
	
		script.close()
		
		i += 1

	print "\nScript has been written.\n"
			
def flagstat(mapping_directory, email, output_file):
	"""
	Using Samtools's flagstat for bam statistics.
	"""
	
	if os.path.isfile(output_file):
		os.remove(output_file)
	
	script = open(output_file,"a")
	script.write("#!/bin/bash\n\n")
	script.close()

	# For PE
	mapping_repo = os.path.join(mapping_directory, "*.bam")

	i = 0
	while i<len(glob.glob(mapping_repo)):

		sample = sorted(glob.glob(mapping_repo))[i]	
		sample_name = sample.split("/")[-1].split(".bam")[0]
		print sample
		
	
		script = open(output_file,"a")
		
		script.write("echo 'module add UHTS/Analysis/samtools/1.4 && samtools flagstat {}' | bsub -u {} -N -J flagstat_{} -o {}\n\n".format(sample, email, sample_name, os.path.join(mapping_directory, sample_name + ".flagstat.txt")))		
	
		script.close()

		
		i += 1


	print "\nScript has been written.\n"
												
def map_filter(mapping_directory, email, output_file):
	"""
	Using Samtools to filter only UNMAPPED (and MATE paired for PE).
	
	https://broadinstitute.github.io/picard/explain-flags.html

	"""
	
	if os.path.isfile(output_file):
		os.remove(output_file)
	
	script = open(output_file,"a")
	script.write("#!/bin/bash\n\n")
	script.close()

	# For PE
	mapping_repo = os.path.join(mapping_directory, "*.bam")
	filter_path = os.path.join(mapping_directory, "filter")
	
	i = 0
	while i<len(glob.glob(mapping_repo)):

		sample = sorted(glob.glob(mapping_repo))[i]	
		file_name = sample.split("/")[-1]
		sample_type = file_name.split(".")[-2]	#all (single) or pair12
		
		if sample_type == "all":	#single 180
		
			sample_name = ".".join((file_name.split(".all")[0].split("."))[0:-1])
			insert_size = (file_name.split(".all")[0].split("."))[-1]
	
			#print sample
			#f: condition valide, F: inverse de condition
			script = open(output_file,"a")
			
			script.write("echo 'mkdir -p {} && module add UHTS/Analysis/samtools/1.3 && samtools view -h -f 4 {} | samtools view -bS - > {}' | bsub -q bgee -m 'cpt133' -M 5000000 -R 'span[ptile=1] rusage[mem=5000]' -u {} -N -J map_filter_{}\n\n".format(filter_path, sample, os.path.join(filter_path, sample_name + ".no_contaminant." + insert_size + ".all.bam"), email, insert_size))		
					
			script.close()
		
		else:
		
			sample_name = ".".join((file_name.split(".pair12")[0].split("."))[0:-1])
			insert_size = (file_name.split(".pair12")[0].split("."))[-1]
	
			#print sample
		
			script = open(output_file,"a")
			
			script.write("echo 'mkdir -p {} && module add UHTS/Analysis/samtools/1.3 && samtools view -h -f 12 {} | samtools view -bS - > {}' | bsub -q bgee -m 'cpt133' -M 5000000 -R 'span[ptile=1] rusage[mem=5000]' -u {} -N -J map_filter_{}\n\n".format(filter_path, sample, os.path.join(filter_path, sample_name + ".no_contaminant." + insert_size + ".pair12.bam"), email, insert_size))		
		
			script.close()
		
		i += 1

	print "\nScript has been written.\n"

def bamtofastq(mapping_directory, email, output_file):
	"""
	Using Picard tools to extract fastq from bam files.
		
	https://broadinstitute.github.io/picard/command-line-overview.html#SamToFastq
	
	"""
	
	if os.path.isfile(output_file):
		os.remove(output_file)
	
	script = open(output_file,"a")
	script.write("#!/bin/bash\n\n")
	script.close()

	# For PE
	mapping_repo = os.path.join(mapping_directory, "*.bam")
	fastq_path = os.path.join(mapping_directory, "fastq")
	
	i = 0
	while i<len(glob.glob(mapping_repo)):

		sample = sorted(glob.glob(mapping_repo))[i]	
		file_name = sample.split("/")[-1]
		sample_type = file_name.split(".")[-2]	#all (single) or pair12
		
		if sample_type == "all":	#single 180
		
			sample_name = ".".join((file_name.split(".all")[0].split("."))[0:-1])
			insert_size = (file_name.split(".all")[0].split("."))[-1]
	
			#print sample
		
			script = open(output_file,"a")
			
			script.write("echo 'mkdir -p {} && module add UHTS/Analysis/picard-tools/2.2.1 && picard-tools SamToFastq I={} F={}' | bsub -q dee-chapuisat -M 20000000 -R 'span[ptile=1] rusage[mem=20000]' -u {} -N -J samtofq_{} -e {}\n\n".format(fastq_path, sample, os.path.join(fastq_path, sample_name + "." + insert_size + ".all.fastq"), email, sample_name, os.path.join(fastq_path, sample_name + ".samtofastq.output.log")))		
		
			script.close()
		
		else:
		
			sample_name = ".".join((file_name.split(".pair12")[0].split("."))[0:-1])
			insert_size = (file_name.split(".pair12")[0].split("."))[-1]
	
			#print sample
		
			script = open(output_file,"a")
			
			script.write("echo 'mkdir -p {} && module add UHTS/Analysis/picard-tools/2.2.1 && picard-tools SamToFastq I={} F={} F2={}' | bsub -q dee-chapuisat -M 20000000 -R 'span[ptile=1] rusage[mem=20000]' -u {} -N -J samtofq_{} -e {}\n\n".format(fastq_path, sample, os.path.join(fastq_path, sample_name + "." + insert_size + ".pair1.fastq"), os.path.join(fastq_path, sample_name + "." + insert_size + ".pair2.fastq"), email, sample_name, os.path.join(fastq_path, sample_name + ".samtofastq.output.log")))	
		
			script.close()
		
		i += 1

	print "\nScript has been written.\n"
	
def fastq_stat(reads_directory, email, output_file):
	"""
	Using BBmap reformat.sh to extract fastq statistics.
	"""
	
	if os.path.isfile(output_file):
		os.remove(output_file)
	
	script = open(output_file,"a")
	script.write("#!/bin/bash\n\n")
	script.close()

	# For PE
	reads_repo = os.path.join(reads_directory, "*.fast*")
	
	i = 0
	while i<len(glob.glob(reads_repo)):

		sample = sorted(glob.glob(reads_repo))[i]	
		sample_name = sample.split("/")[-1].split(".fastq.gz")[0]
		#sample_name = sample.split("/")[-1].split(".fasta")[0]
		#print sample
		
	
		script = open(output_file,"a")
		
		script.write("echo 'module add UHTS/Analysis/BBMap/37.82 && reformat.sh in={}' | bsub -u {} -N -J reformat_{} -e {}\n\n".format(sample, email, sample_name, sample_name + ".stat.output.txt"))		
	
		script.close()
		
		i += 1

	print "\nScript has been written.\n"

def fastq_stat_table(fastq_stat_directory, output_file):

	"""
	Statistics after fastq_stat
	"""
	
	if os.path.isfile(output_file):
		os.remove(output_file)
	
	"""
	script = open("1_Tdi_is_3000.sorted.filtered.pair1.stat.output.txt","r")
	stat = [line.split() for line in script if re.search("bases",line)][0][3]
	print stat
	script.close()

	"""	
	data_repo = os.path.join(fastq_stat_directory, "*.stat.output.txt")
	
	script = open(output_file,"a")
	script.write("Sample\tReads number\tBases number\n")
	sum_base = 0
	
	i = 0
	while i<len(glob.glob(data_repo)):
	
		data = sorted(glob.glob(data_repo))[i]
		sample_name = (data.split("/")[-1]).split(".")[0]
			
		fastq_stat_output_log = open(data,"r")
		
		stat = [line.split() for line in fastq_stat_output_log if re.search("bases",line)][0]
		sum_base += int(stat[3])
		
		fastq_stat_output_log.close()
		
		script.write("{}\t{}\t{}\n".format(sample_name, stat[1], stat[3]))
				
		i += 1
	
	script.close()

	print "Table has been written."
	print "Sum bases = ", sum_base
		
def flag_fastq_stat(bwa_directory, reads_directory, output_file):

	"""
	Statistics for Paired-end after trimmomatic().
	"""
	
	if os.path.isfile(output_file):
		os.remove(output_file)
		
	bwa_repo = os.path.join(bwa_directory, "*.flagstat.txt")
	
	script = open(output_file,"a")
	script.write("Lib\tSample\tInput reads (Tot)\tProperly paired mapped (Tot reads)\tFiltered total reads\tFiltered total base\tCoverage\n")
	
	i = 0
	while i<len(glob.glob(bwa_repo)):
	
		flagstat_file = sorted(glob.glob(bwa_repo))[i]	
		#print flagstat_file #bwa_eukaryote/Sm2.nxtrim.3000.pair12.flagstat.txt
		flagstat_decompose = flagstat_file.split("/")[-1].split(".")
		
		fastq_stat_file_p1 = os.path.join(reads_directory, flagstat_decompose[0] + "." + flagstat_decompose[1].replace("trim", "trim.filtered") + "." + flagstat_decompose[2] + ".pair1.stat.error.txt")
		
		fastq_stat_file_p2 = os.path.join(reads_directory, flagstat_decompose[0] + "." + flagstat_decompose[1].replace("trim", "trim.filtered") + "." + flagstat_decompose[2] + ".pair2.stat.error.txt")
							
		#print flagstat_file	
		#print fastq_stat_file_p1
		#print fastq_stat_file_p2

		flagstat_output_log = open(flagstat_file,"r")
		input_read = ([line.split() for line in flagstat_output_log if re.search("paired in sequencing",line)][0])
		flagstat_output_log.close()
		
		flagstat_output_log = open(flagstat_file,"r")
		properly_paired = ([line.split() for line in flagstat_output_log if re.search("properly paired",line)][0])
		flagstat_output_log.close()
		
		fastq_stat_p1_output_log = open(fastq_stat_file_p1,"r")
		stat_p1 = ([line.split() for line in fastq_stat_p1_output_log if re.search("Input:",line)][0])
		fastq_stat_p1_output_log.close()

		fastq_stat_p2_output_log = open(fastq_stat_file_p2,"r")
		stat_p2 = ([line.split() for line in fastq_stat_p2_output_log if re.search("Input:",line)][0])
		fastq_stat_p2_output_log.close()
				
		
		script.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(flagstat_decompose[2], flagstat_decompose[0], input_read[0], properly_paired[0] + " " + properly_paired[5] + ")", int(stat_p1[1]) + int(stat_p2[1]), int(stat_p1[3]) + int(stat_p2[3]), (int(stat_p1[3]) + int(stat_p2[3])) / 300000000))
		
		i += 1
	
	script.close()

	print "Table has been written."
	
def mark_duplicates(mapping_directory, email, output_file):
	"""
	Using Picard tools to remove duplicates (PCR and optical).
	https://broadinstitute.github.io/picard/command-line-overview.html	
	"""
	
	if os.path.isfile(output_file):
		os.remove(output_file)
	
	script = open(output_file,"a")
	script.write("#!/bin/bash\n\n")
	script.close()

	# For PE
	mapping_repo = os.path.join(mapping_directory, "*.bam")
	duplicate_path = os.path.join(mapping_directory, "mark_duplicates")
	
	i = 0
	while i<len(glob.glob(mapping_repo)):

		sample = sorted(glob.glob(mapping_repo))[i]	
		file_name = sample.split("/")[-1]
		sample_name = ".".join((file_name.split(".pair12")[0].split("."))[0:-1])
		insert_size = (file_name.split(".pair12")[0].split("."))[-1]

		#print sample
			
		script = open(output_file,"a")
		
		script.write("echo 'mkdir -p {} && module add UHTS/Analysis/picard-tools/2.2.1 && picard-tools MarkDuplicates I={} O={} REMOVE_DUPLICATES=true M={}' | bsub -q dee-hugemem -M 50000000 -R 'span[ptile=1] rusage[mem=50000]' -u {} -N -J mark_duplicates_{} -e {}\n\n".format(duplicate_path, sample, os.path.join(duplicate_path, sample_name + ".nodup." + insert_size + ".pair12.bam"), os.path.join(duplicate_path, sample_name + ".nodup." + insert_size + ".stat.txt"), email, sample_name, os.path.join(duplicate_path, sample_name + ".mark_duplicates.log")))			
		script.close()
		
		i += 1

	print "\nScript has been written.\n"

def rename(fasta_file, output_file):
	"""
	Rename header:
	Initial: >scaffold24|size30367
	After: >scaffold24_size30367
	"""
	
	if os.path.isfile(output_file):
		os.remove(output_file)
	
	seq_record_ass = []
	i = 1
	
	for seq_record in SeqIO.parse(fasta_file, "fasta"):
		split_id = seq_record.id.split("|")
		#seq_record.id = "Gap" + str(i) + "_" + str(len(seq_record))
		seq_record.id = "Contig" + str(i) + "_" + str(len(seq_record))
		
			
		seq_record_ass.append(seq_record)
		
		i += 1
	
	SeqIO.write(seq_record_ass, output_file, "fasta")

def rename_TE(fasta_file, output_file):
	"""
	Rename header:
	Initial: >DNA_comp4395_g1_i1#DNA/DNA
	After: >DNA_comp4395_g1_i1
	"""
	
	if os.path.isfile(output_file):
		os.remove(output_file)
	
	seq_record_ass = []
	i = 1
	
	for seq_record in SeqIO.parse(fasta_file, "fasta"):
		split_id = seq_record.id.split("#")
		seq_record.id = split_id[0]
		seq_record.description = ""
		
			
		seq_record_ass.append(seq_record)
		
		i += 1
	
	SeqIO.write(seq_record_ass, output_file, "fasta")
	
def reformat(fasta_file, output_file):
	"""
	Reformat with biopython
	"""
	
	if os.path.isfile(output_file):
		os.remove(output_file)
	
	seq_record_ass = []
	
	for seq_record in SeqIO.parse(fasta_file, "fasta"):
			
		seq_record_ass.append(seq_record)
	
	SeqIO.write(seq_record_ass, output_file, "fasta")
	
def busco_assignment(busco_tsv, tax_list, output):
	"""
	Busco assignment
	"""

	if os.path.isfile(output):
		os.remove(output)
				
	busco_table_file = open(busco_tsv,"r")
	for busco in busco_table_file:
	
		if '#' in busco:
			pass
		else:
			info = re.split("[\r\t\n]",busco)
			#print info
			
			if info[1] != "Missing":

				phylums = open(tax_list,"r")
				for phylum in phylums:
				
					info_phylum = re.split("[\r\t\n]",phylum) 
					
					if info_phylum[0] == info[2]:
					
						script = open(output,"a")			
						script.write("{}\t{}\t{}\t{}\n".format(info[0], info[1], info[2], info_phylum[1]))
						script.close()
						
						break

				phylums.close()
			
			
	busco_table_file.close()

	print "\n\nTable has been written.\n"

def cegma_assignment(cegma_gff, tax_list, output):
	"""
	Cegma assignment
	"""

	if os.path.isfile(output):
		os.remove(output)
	
	list_scaffold = []
				
	cegma_table_file = open(cegma_gff,"r")
	for cegma in cegma_table_file:
	
		if '#' in cegma:
			pass
		else:
			info = re.split("[\r\t\n]",cegma)
			
			if info[0] in list_scaffold:
				pass
			else:
			
				list_scaffold.append(info[0])
			
				phylums = open(tax_list,"r")
				for phylum in phylums:
				
					info_phylum = re.split("[\r\t\n]",phylum) 
					
					if info_phylum[0] == info[0]:
					
						script = open(output,"a")			
						script.write("{}\t{}\n".format(info[0], info_phylum[1]))
						script.close()
						
						break

			phylums.close()
			
			
	cegma_table_file.close()

	print "\n\nTable has been written.\n"
	
def busco_coverage(busco_tsv, coverage_file, fasta_file, output):
	"""
	Busco coverage
	"""

	if os.path.isfile(output):
		os.remove(output)
	
	i = 0
				
	buscos = open(busco_tsv,"r")
	
	for busco in buscos:
	
		#print re.split("[\r\n]",busco)
		
		if '#' in busco:
			pass
		else:
			#info = re.split("[\r\t\n]",busco)
			#print info
			busco_status = re.split("[\r\t\n]",busco)[1]
			busco_scaffold = re.split("[\r\t\n]",busco)[2]
			
		 	print busco_scaffold
		 	
			if busco_status != "Missing":

				for seq_record in SeqIO.parse(fasta_file, "fasta"):
				
					if busco_scaffold == seq_record.id:	# same scaffolds name
					
						length_scaffold = len(seq_record)						
						break
								
				coverages = open(coverage_file,"r")
				
				for coverage in coverages:
				
					coverage_scaffold = re.split("[\r\t\n]",coverage)[0]
					coverage_nb = re.split("[\r\t\n]",coverage)[1]
					
					if busco_scaffold == coverage_scaffold:
						break
				
				coverages.close()
				
				
				if busco_status == "Duplicated":
					color = 1
				if busco_status == "Complete":
					color = 2		
				if busco_status == "Fragmented":
					color = 3
			
				script = open(output,"a")			
				script.write("{}\t{}\t{}\t{}\n".format(re.split("[\r\n]",busco)[0], length_scaffold, coverage_nb, color))
				script.close()
						
				i += 1
				
				print i
				
	buscos.close()
					
	print "\n\nTable has been written.\n"

def cegma_reformat(cegma_id_file, cegma_gff_file, output):
	"""
	cegma_reformat
	"""

	if os.path.isfile(output):
		os.remove(output)
	
	i = 0
				
	cegma_ids = open(cegma_id_file,"r")
	
	for cegma_id in cegma_ids:
	
		#print re.split("[\r\n]",busco		
		
		ko_name_id = re.split("[\r\t\n]",cegma_id)[0]
	 	
		cegma_gffs = open(cegma_gff_file,"r")
		
		for cegma_gff in cegma_gffs:
		
			ko_name_gff = re.split("[\r\t\n]",cegma_gff)[8]
			scaffold_gff = re.split("[\r\t\n]",cegma_gff)[0]
		
			if ko_name_id == ko_name_gff:
				break
		
		cegma_gffs.close()
		
		script = open(output,"a")			
		script.write("{}\t{}\n".format(ko_name_id, scaffold_gff))
		script.close()
				

	cegma_ids.close()					
						
	print "\n\nTable has been written.\n"

def cegma_coverage(cegma_summary, coverage_file, fasta_file, output):
	"""
	Cegma coverage
	"""

	if os.path.isfile(output):
		os.remove(output)
	
	i = 0
				
	cegmas = open(cegma_summary,"r")
	
	for cegma in cegmas:
	
		#print re.split("[\r\n]",busco)
		
		cegma_scaffold = re.split("[\r\t\n]",cegma)[1]
		
	 	print cegma_scaffold
	 	
		for seq_record in SeqIO.parse(fasta_file, "fasta"):
		
			if cegma_scaffold == seq_record.id:	# same scaffolds name
			
				length_scaffold = len(seq_record)						
				break
						
		coverages = open(coverage_file,"r")
		
		for coverage in coverages:
		
			coverage_scaffold = re.split("[\r\t\n]",coverage)[0]
			coverage_nb = re.split("[\r\t\n]",coverage)[1]
			
			if cegma_scaffold == coverage_scaffold:
				break
		
		coverages.close()
		
		script = open(output,"a")			
		script.write("{}\t{}\t{}\t{}\n".format(re.split("[\r\n]",cegma)[0], length_scaffold, coverage_nb, 4))
		script.close()
				
		i += 1
		
		print i
		
	cegmas.close()
					
	print "\n\nTable has been written.\n"
	
	
	
def diff_scaf_busco_assign(scaf_file, busco_assignemnt_file):
	"""
	Scaffold_list after blobtool_contaminants
	"""
	
	buscos = open(busco_assignemnt_file,"r")
		
	for busco in buscos:
		
		scaffolds = open(scaf_file,"r")
		
		exist = False
		
		for scaffold in scaffolds:
		
			if re.split("[\r\t\n]",busco)[2] == re.split("[\r\t\n]",scaffold)[0]:				
				exist = True
				break
					
		scaffolds.close()
		
		if exist == False:
		
			print busco
			
	buscos.close()	


def busco_cegma_cov(busco_coverage_file, cegma_coverage_file, output_file):
	"""
	Coverage with color depending on type for R
	"""
	
	busco_dic = {}
	cegma_dic = {}
	
	busco_list = []
	cegma_in_busco = []
	
	buscos = open(busco_coverage_file,"r")
		
	for busco in buscos:
		
		busco_scaf_name = re.split("[\r\t\n]",busco)[2]
		
		if busco_scaf_name in busco_dic:
			busco_dic[busco_scaf_name] += 1
		else: 
			busco_dic[busco_scaf_name] = 1			

	buscos.close()
	
	cegmas = open(cegma_coverage_file,"r")
		
	for cegma in cegmas:
		
		cegma_scaf_name = re.split("[\r\t\n]",cegma)[1]
		
		if cegma_scaf_name in cegma_dic:
			cegma_dic[cegma_scaf_name] += 1
		else: 
			cegma_dic[cegma_scaf_name] = 1
					
	cegmas.close()
	
	
	buscos = open(busco_coverage_file,"r")
	
	for busco in buscos:
		
		busco_scaf_name = re.split("[\r\t\n]",busco)[2]
		
		if busco_dic[busco_scaf_name] == 1 :	# unique busco id
		
			if busco_scaf_name in cegma_dic:	#shared with cegma dic
			
				if cegma_dic[busco_scaf_name] == 1:	# 1 busco id is shared with 1 cegma id
				
					script = open(output_file,"a")			
					script.write("{}\t{}\t{}\t{}\t{}\n".format(busco_scaf_name, int(re.split("[\r\t\n]",busco)[7]), float(re.split("[\r\t\n]",busco)[8]), str(int(re.split("[\r\t\n]",busco)[9]) + 4), "Dot"))
					script.close()
					
					cegma_in_busco.append(busco_scaf_name)
					
				else:	# 1 busco id is shared with >1 cegma id
				
					script = open(output_file,"a")			
					script.write("{}\t{}\t{}\t{}\t{}\n".format(busco_scaf_name, int(re.split("[\r\t\n]",busco)[7]), float(re.split("[\r\t\n]",busco)[8]), "9", "Cross"))
					script.close()
				
					cegma_in_busco.append(busco_scaf_name)
					
			else:	# 1 busco id with 0 cegma id
			
				script = open(output_file,"a")			
				script.write("{}\t{}\t{}\t{}\t{}\n".format(busco_scaf_name, int(re.split("[\r\t\n]",busco)[7]), float(re.split("[\r\t\n]",busco)[8]), str(re.split("[\r\t\n]",busco)[9]), "Dot"))
				script.close()
				
		else:	# n busco id in the same scaffold
		
			if busco_scaf_name in cegma_dic:	# n busco id and n cegma id

				script = open(output_file,"a")			
				script.write("{}\t{}\t{}\t{}\t{}\n".format(busco_scaf_name, int(re.split("[\r\t\n]",busco)[7]), float(re.split("[\r\t\n]",busco)[8]), "9", "Cross"))
				script.close()
				
				cegma_in_busco.append(busco_scaf_name)
				
			else:	# n busco id and 0 cegma id
			
				script = open(output_file,"a")			
				script.write("{}\t{}\t{}\t{}\t{}\n".format(busco_scaf_name, int(re.split("[\r\t\n]",busco)[7]), float(re.split("[\r\t\n]",busco)[8]), "8", "Cross"))
				script.close()
			
	buscos.close()


	cegmas = open(cegma_coverage_file,"r")
	
	for cegma in cegmas:
		
		cegma_scaf_name = re.split("[\r\t\n]",cegma)[1]
		
		if cegma_scaf_name not in cegma_in_busco:	# cegma that are not shared with busco
		
			if cegma_dic[cegma_scaf_name] == 1:	# 1 cegma id and 0 busco id
			
				script = open(output_file,"a")			
				script.write("{}\t{}\t{}\t{}\t{}\n".format(cegma_scaf_name, int(re.split("[\r\t\n]",cegma)[2]), float(re.split("[\r\t\n]",cegma)[3]), str(re.split("[\r\t\n]",cegma)[4]), "Dot"))
				script.close()
				
			else: #n cegma id and 0 busco id
			
				script = open(output_file,"a")			
				script.write("{}\t{}\t{}\t{}\t{}\n".format(cegma_scaf_name, int(re.split("[\r\t\n]",cegma)[2]), float(re.split("[\r\t\n]",cegma)[3]), str(re.split("[\r\t\n]",cegma)[4]), "Cross"))
				script.close()
				
								
	cegmas.close()
			
	

def random_coverage(fasta_file, output_file):
	"""
	Random coverage for blobtools
	"""
	
	if os.path.isfile(output_file):
		os.remove(output_file)
			

	
	for seq_record in SeqIO.parse(fasta_file, "fasta"):
		
		script = open(output_file,"a")			
		script.write("{}\t{}\n".format(seq_record.id, random.randint(1, 56)))
		script.close()
				

def scaffolds_busco_cegma(busco_summary, cegma_summary, contaminant_scaffolds, output):
	"""
	Extract list of scaffolds that have busco and/or cegma id. 
	Then you have to use bbmap to extract scaffolds sequence
	
	busco_summary: full_table_output.tsv
	cegma_summary: output.cegma.gff
	contaminant_scaffolds: scaf_contaminant.txt after blobtool_contaminants()
	"""

	if os.path.isfile(output):
		os.remove(output)
	
	busco_ids = []
	cegma_ids = []
	contaminant_ids = []
	busco_cegma_ids = []
				
	busco_table_file = open(busco_summary,"r")
	
	for busco in busco_table_file:
	
		if '#' in busco:
			pass
		else:
			info_busco = re.split("[\r\t\n]",busco)
			#print info
			
			if info_busco[1] != "Missing":

				busco_ids.append(info_busco[2])
				busco_cegma_ids.append(info_busco[2])
						
	busco_table_file.close()


	cegma_table_file = open(cegma_summary,"r")
	
	for cegma in cegma_table_file:

		info_cegma = re.split("[\r\t\n]",cegma)	 	
		cegma_ids.append(info_cegma[0])
		busco_cegma_ids.append(info_cegma[0])

	cegma_table_file.close()	
	
	
	contaminants_table_file = open(contaminant_scaffolds,"r")

	for contaminants in contaminants_table_file:

		info_contaminants = re.split("[\r\t\n]",contaminants)	 	
		contaminant_ids.append(info_contaminants[0])

	contaminants_table_file.close()	
	
	
	ids_not_cont = 0
	
	script = open(output,"a")
	
	for ids in list(set(busco_cegma_ids)):
	
		if ids not in contaminant_ids:
		
			ids_not_cont += 1
			script.write("{}\n".format(ids))
	
	script.close()
	
	print "Unique Busco: ", len(list(set(busco_ids)))
	print "Unique Cegma: ", len(list(set(cegma_ids)))
	print "Unique Busco and Cegma: ", len(list(set(busco_cegma_ids)))
	print "Unique Busco and Cegma not contaminated: ", ids_not_cont
						
	print "\n\nTable has been written.\n"
				

def cat_reseq(raw_read_directory, email, output_file):

	"""
	raw_read_directory = /archive/dee/schwander/ptranvan/mites/resequence/350 (with the insert size)
	cat reseq
	"""
		
	if os.path.isfile(output_file):
		os.remove(output_file)
		
	script = open(output_file,"a")
	script.write("#!/bin/bash\n\n")
	script.close()
	
	"""
	dic_350 = {"128": "On1", "129": "On3", "130": "On4", "131": "On5", "132": "On7", "133": "Os1", "134": "Os2", "135": "Os3", "136": "Os4", "137": "Os7", "138": "As1", "139": "As2", "140": "As4", "141": "As5", "142": "As7", "143": "Sm1", "144": "Sm3", "145": "Sm4", "146": "Sm5", "147": "Sm6", "148": "Np1", "149": "Np2", "150": "Np3"}
	
	for key, value in sorted(dic_350.iteritems(), key=lambda (k,v): (v,k)):
	
		print value
	
		script = open(output_file,"a")
				
		script.write("echo 'mkdir -p 350 && zcat {} > {}' | bsub -u {} -N -J zcat_{}_R1 \n\n".format(os.path.join(raw_read_directory, "*{}*R1*".format(key)), os.path.join("350", "{}.raw.350.pair1.fastq".format(value)), email, key))

		script.write("echo 'mkdir -p 350 && zcat {} > {}' | bsub -u {} -N -J zcat_{}_R2 \n\n".format(os.path.join(raw_read_directory, "*{}*R2*".format(key)), os.path.join("350", "{}.raw.350.pair2.fastq".format(value)), email, key))
								
		script.close()
	"""
	
	dic_550 = {"151": "On1", "152": "On3", "153": "On4", "154": "On5", "155": "On7", "156": "Os1", "157": "Os2", "158": "Os3", "159": "Os4", "160": "Os7", "161": "As1", "162": "As2", "163": "As4", "164": "As5", "165": "As7", "166": "Sm1", "167": "Sm3", "168": "Sm4", "169": "Sm5", "170": "Sm6", "171": "Np1", "172": "Np2", "173": "Np3"}
	
	for key, value in sorted(dic_550.iteritems(), key=lambda (k,v): (v,k)):
	
		print value
	
		script = open(output_file,"a")
				
		script.write("echo 'mkdir -p 550 && zcat {} > {}' | bsub -u {} -N -J zcat_{}_R1 \n\n".format(os.path.join(raw_read_directory, "*{}*R1*".format(key)), os.path.join("550", "{}.raw.550.pair1.fastq".format(value)), email, key))

		script.write("echo 'mkdir -p 550 && zcat {} > {}' | bsub -u {} -N -J zcat_{}_R2 \n\n".format(os.path.join(raw_read_directory, "*{}*R2*".format(key)), os.path.join("550", "{}.raw.550.pair2.fastq".format(value)), email, key))
								
		script.close()
						
	
	print "done"
		#print "%s: %s" % (key, value)

def bwa_reseq(trim_directory, db, specie, output_file):
	"""
	Map reads 
	"""
	
	output_dir = "bwa"
	
	if os.path.isfile(output_file):
		os.remove(output_file)
	
	script = open(output_file,"a")
	script.write("#!/bin/bash\n\n")
	script.write("module add UHTS/Aligner/bwa/0.7.13\n\n")
	script.write("mkdir bwa\n\n")
	script.close()


	data_repo = os.path.join(trim_directory, "*")
	
	i = 0
	while i<len(glob.glob(data_repo)):
	
		ins_length_path = sorted(glob.glob(data_repo))[i]	
		read_repo = os.path.join(ins_length_path, "{}*.pair1.fastq",format(specie))
		
		j = 0
		while j<len(glob.glob(read_repo)):
		
			trim_data = sorted(glob.glob(read_repo))[j]
			#print pe_data # trim/350/160603_SND393_A_L002_HYI-120_output_log.txt
			
			insert_size = trim_data.split("/")[-2]
			#print insert_size
						
			sample_name = (trim_data.split("/")[-1]).split(".")[0]
			#print sample_name
			
			script = open(output_file,"a")	
								
			script.write("bwa mem -M {} {} {} -t 15 > {} 2> {}\n\n".format(output_dir, db, trim_data, trim_data.replace("pair1", "pair2"), os.path.join(output_dir, trim_data.split("pair1")[0] + "pair12.sam"), os.path.join(output_dir, trim_data.split("pair1")[0] + "pair12.bwa.output.log")))	
				
			script.close()	
	
			j+=1
			
		i += 1
	
	script.close()
	
	print "\nScript has been written.\n"

def coverage_manual(stat_parse, contaminant_scaffolds, output_file):
	"""
	stat_parse.txt
	contaminant_scaffolds.txt 
	"""
		
	if os.path.isfile(output_file):
		os.remove(output_file)
	
	contaminant_list = []
	
	contaminant_scaffolds_file = open(contaminant_scaffolds,"r")
	
	for contaminant in contaminant_scaffolds_file:
	
		contaminant_list.append(re.split("[\r\t\n]",contaminant)[0])
						
	contaminant_scaffolds_file.close()

	#print contaminant_list
	
	stat_parse_file = open(stat_parse,"r")
	
	for coverage in stat_parse_file:
	
		scaffold_name = re.split("[\r\t\n]",coverage)[0]
		
		if scaffold_name in contaminant_list:
			pass
		else:
			
			script = open(output_file,"a")									
			script.write("{}".format(coverage))				
			script.close()	

	stat_parse_file.close()

	
	print "\nScript has been written.\n"

def blast_manual(blast, contaminant_scaffolds, output_file):
	"""
	sp.vs.nt.cul5.1e25.megablast.out
	contaminant_scaffolds.txt 
	"""
		
	if os.path.isfile(output_file):
		os.remove(output_file)
	
	contaminant_list = []
	
	contaminant_scaffolds_file = open(contaminant_scaffolds,"r")
	
	for contaminant in contaminant_scaffolds_file:
	
		contaminant_list.append(re.split("[\r\t\n]",contaminant)[0])
						
	contaminant_scaffolds_file.close()

	#print contaminant_list
	
	blast_file = open(blast,"r")
	
	for coverage in blast_file:
	
		scaffold_name = re.split("[\r\t\n]",coverage)[0]
		
		if scaffold_name in contaminant_list:
			pass
		else:
			
			script = open(output_file,"a")									
			script.write("{}".format(coverage))				
			script.close()	

	blast_file.close()

	
	print "\nScript has been written.\n"

def read_size(source_directory, output_file):

	"""
	>ls source_directory
	1_Tdi  1_Tps  2_Tcm 
	>ls 1_Tdi
	fastq/
		
	"""
		
	if os.path.isfile(output_file):
		os.remove(output_file)
	
	script = open(output_file,"a")
	script.write("#!/bin/bash\n\n")
	script.close()
	

	data_repo = os.path.join(source_directory, "*")
	
	commands = "awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}'"
	i = 0
	while i<len(glob.glob(data_repo)):
	
		specie_path = sorted(glob.glob(data_repo))[i]	
		specie = specie_path.split("/")[-1]
					
		read_repo = os.path.join(os.path.join(specie_path, "fastq"), "*fastq")
		

		j = 0
		while j<len(glob.glob(read_repo)):
		
		
			pe_data = sorted(glob.glob(read_repo))[j]
			sample_name = (pe_data.split("/")[-1]).split(".")[0]
		
			script = open(output_file,"a")
			
			script.write("{} {} > {}.size_distribution.txt\n\n".format(commands, pe_data, os.path.join("/".join(pe_data.split("/")[0:-1]), sample_name)))
					
			script.close()

			j+=1

		i += 1
		
				
	print "Script has been written."
											
																								
def main(argv):
	
	mod=[]
	
	mod.append('\n%(prog)s -s fastuniq -i1 <raw_read_merged_directory> -o <output_file>')
	mod.append('%(prog)s -s nxtrim -i1 <raw_read_directory> -e <email> -o <output_file>')
	mod.append('%(prog)s -s cat_nextclip -i1 <nextclip_directory> -e <email> -o <output_file>')
	mod.append('%(prog)s -s trimmomatic -i1 <raw_read_directory> -i2 <nextclip_directory> -e <email> -o <output_file>')
	mod.append('%(prog)s -s MP_stat -i1 <nxtrim_directory> -o <output_file>')
	mod.append('%(prog)s -s PE_stat -i1 <trim_directory> -o <output_file>')
	mod.append('%(prog)s -s blobtool_contaminants -i1 <blob_table> -e <email> -o <output_file>')
	mod.append('%(prog)s -s bowtie2 -i1 <reads.list> -i2 <db> -e <email> -o <output_file>')
	mod.append('%(prog)s -s bwa -i1 <reads.list> -i2 <db> -e <email> -o <output_file>')
	mod.append('%(prog)s -s bbmap -i1 <reads.list> -e <email> -o <output_file>')
	mod.append('%(prog)s -s bbmap_sam -i1 <reads.list> -e <email> -o <output_file>')
	mod.append('%(prog)s -s samtofastq -i1 <mapping_directory> -e <email> -o <output_file>')
	mod.append('%(prog)s -s sort -i1 <mapping_directory> -e <email> -o <output_file>')
	mod.append('%(prog)s -s sort_query -i1 <mapping_directory> -e <email> -o <output_file>')
	mod.append('%(prog)s -s flagstat -i1 <mapping_directory> -e <email> -o <output_file>')
	mod.append('%(prog)s -s map_filter -i1 <mapping_directory> -e <email> -o <output_file>')
	mod.append('%(prog)s -s bamtofastq -i1 <mapping_directory> -e <email> -o <output_file>')
	mod.append('%(prog)s -s fastq_stat -i1 <read_directory> -e <email> -o <output_file>')
	mod.append('%(prog)s -s fastq_stat_table -i1 <fastq_stat_directory> -o <output_file>')
	mod.append('%(prog)s -s flag_fastq_stat -i1 <bwa_directory> -i2 <read_directory> -o <output_file>')
	mod.append('%(prog)s -s mark_duplicates -i1 <mapping_directory> -e <email> -o <output_file>')
	mod.append('%(prog)s -s rename -i1 <fasta_file> -o <output_file>')
	mod.append('%(prog)s -s rename_TE -i1 <fasta_file> -o <output_file>')
	mod.append('%(prog)s -s busco_assignment -i1 <busco_tsv> -i2 <phylum_pre> -o <output_file>')
	mod.append('%(prog)s -s busco_coverage -i1 <busco_tsv> -i2 <coverage_file> -i3 <fasta_file> -o <output_file>')
	mod.append('%(prog)s -s cegma_reformat -i1 <cegma_id_file> -i2 <cegma_gff_file> -o <output_file>')
	mod.append('%(prog)s -s cegma_coverage -i1 <cegma_summary> -i2 <coverage_file> -i3 <fasta_file> -o <output_file>')
	mod.append('%(prog)s -s busco_cegma_cov -i1 <busco_coverage> -i2 <cegma_coverage> -o <output_file>')
	mod.append('%(prog)s -s random_coverage -i1 <fasta_file> -o <output_file>')
	mod.append('%(prog)s -s scaffolds_busco_cegma -i1 <busco_summary> -i2 <cegma_summary> -i3 <contaminant_scaffolds> -o <output_file>')
	mod.append('%(prog)s -s bwa_reseq -i1 <trim_directory> -i2 <db> -i3 <specie> -o <output_file>')
	mod.append('%(prog)s -s coverage_manual -i1 <stat_parse> -i2 <contaminant_scaffolds> -o <output_file>')
	mod.append('%(prog)s -s read_size -i1 <source_directory> -o <output_file>')


	parser = argparse.ArgumentParser(prog = 'asv3.py',
                                 usage = "\n".join(mod))

	parser.add_argument('-s', action='store', dest='step_value',
	                    help='Step')
	                                                 	
	parser.add_argument('-i1', action='store', dest='input_value',
	                    help='Input 1')

	parser.add_argument('-i2', action='store', dest='input2_value',
	                    help='Input 2')

	parser.add_argument('-i3', action='store', dest='input3_value',
	                    help='Input 3')
	                    
	parser.add_argument('-e', action='store', dest='email_value',
	                    help='Email')	                
	                    	                    	
	parser.add_argument('-o', action='store', dest='output_value',
	                    help='Output')
	                    	                    	
	parser.add_argument('--version', action='version', version='%(prog)s 1.0')

	results = parser.parse_args()


	# Run FastUniq
	
	if results.step_value == "fastuniq" and results.input_value and results.output_value:
		fastuniq(results.input_value, results.output_value)

	# Run nextclip in Vital-IT cluster.
	
	if results.step_value == "nxtrim" and results.input_value  and results.email_value and results.output_value:
		nxtrim(results.input_value, results.email_value, results.output_value)

	# Fetch real Mate pair.
	
	if results.step_value == "cat_nxtrim" and results.input_value  and results.email_value and results.output_value:
		cat_nxtrim(results.input_value, results.email_value, results.output_value)

	# Run trimmomatic in Vital-IT cluster.		
				
	if results.step_value == "trimmomatic" and results.input_value and results.input2_value and results.email_value and results.output_value:
		trimmomatic(results.input_value, results.input2_value, results.email_value, results.output_value)

		
	# Create a table for MP's stats.			
				
	if results.step_value == "MP_stat" and results.input_value and results.output_value:
		MP_stat(results.input_value, results.output_value)
				
	# Create a table of PE's stats.			
				
	if results.step_value == "PE_stat" and results.input_value and results.output_value:
		PE_stat(results.input_value, results.output_value)

	# Create a list of scaffolds with potential contaminants.			

	if results.step_value == "taxonomy_list" and results.input_value:
		taxonomy_list(results.input_value)
				
	if results.step_value == "blobtool_contaminants" and results.input_value:
		blobtool_contaminants(results.input_value)

	if results.step_value == "link_contaminants" and results.input_value and results.input2_value:
		link_contaminants(results.input_value, results.input2_value)

	if results.step_value == "contaminants_filtration" and results.input_value and results.output_value:
		contaminants_filtration(results.input_value, results.output_value)

	# Using bwa to map reads on scaffolds.			
				
	if results.step_value == "bwa" and results.input_value and results.input2_value and results.email_value and results.output_value:
		bwa(results.input_value, results.input2_value, results.email_value, results.output_value)

	# Using bbmap to map reads on scaffolds.			
				
	if results.step_value == "bbmap" and results.input_value and results.email_value and results.output_value:
		bbmap(results.input_value, results.email_value, results.output_value)

	if results.step_value == "bbmap_sam" and results.input_value and results.email_value and results.output_value:
		bbmap_sam(results.input_value, results.email_value, results.output_value)
				
	# Using samtools to sort bam files.
				
	if results.step_value == "sort" and results.input_value and results.email_value and results.output_value:
		sort(results.input_value, results.email_value, results.output_value)
				
	if results.step_value == "sort_query" and results.input_value and results.email_value and results.output_value:
		sort_query(results.input_value, results.email_value, results.output_value)
				
	# Using samtools for bam statistics.
				
	if results.step_value == "flagstat" and results.input_value and results.email_value and results.output_value:
		flagstat(results.input_value, results.email_value, results.output_value)
																						
	# Using samtools for bam filtering.
				
	if results.step_value == "map_filter" and results.input_value and results.email_value and results.output_value:
		map_filter(results.input_value, results.email_value, results.output_value)
																										
	# Using Picard tools to extract fastq from bam.
				
	if results.step_value == "bamtofastq" and results.input_value and results.email_value and results.output_value:
		bamtofastq(results.input_value, results.email_value, results.output_value)

	# Using BBmap reformat.sh to extract fastq statistics.
				
	if results.step_value == "fastq_stat" and results.input_value and results.email_value and results.output_value:
		fastq_stat(results.input_value, results.email_value, results.output_value)

	# Table with summary of fastq_stat result
				
	if results.step_value == "fastq_stat_table" and results.input_value and results.output_value:
		fastq_stat_table(results.input_value, results.output_value)

	# Table to check statisttic results from flagstat and reformat stat			
				
	if results.step_value == "flag_fastq_stat" and results.input_value and results.input2_value and results.output_value:
		flag_fastq_stat(results.input_value, results.input2_value, results.output_value)

	# Using Picard tools to remove duplicates from bam files.
				
	if results.step_value == "mark_duplicates" and results.input_value and results.email_value and results.output_value:
		mark_duplicates(results.input_value, results.email_value, results.output_value)

	# Rename header of FASTA.
				
	if results.step_value == "rename" and results.input_value and results.output_value:
		rename(results.input_value, results.output_value)

	# Rename header of FASTA. (TE jens)
				
	if results.step_value == "rename_TE" and results.input_value and results.output_value:
		rename_TE(results.input_value, results.output_value)
		
	# Reformat FASTA.
				
	if results.step_value == "reformat" and results.input_value and results.output_value:
		reformat(results.input_value, results.output_value)
		
	# Create a table for busco assignment.			
				
	if results.step_value == "busco_assignment" and results.input_value and results.input2_value and results.output_value:
		busco_assignment(results.input_value, results.input2_value, results.output_value)

	if results.step_value == "cegma_assignment" and results.input_value and results.input2_value and results.output_value:
		cegma_assignment(results.input_value, results.input2_value, results.output_value)
		
	if results.step_value == "busco_coverage" and results.input_value and results.input2_value and results.output_value:
		busco_coverage(results.input_value, results.input2_value, results.input3_value, results.output_value)

	if results.step_value == "cegma_reformat" and results.input_value and results.input2_value and results.output_value:
		cegma_reformat(results.input_value, results.input2_value, results.output_value)

	if results.step_value == "cegma_coverage" and results.input_value and results.input2_value and results.output_value:
		cegma_coverage(results.input_value, results.input2_value, results.input3_value, results.output_value)
		
			
	# Extract fasta from scaffolds after blobtools_contaminant()			

			
	if results.step_value == "diff_scaf_busco_assign" and results.input_value and results.input2_value:
		diff_scaf_busco_assign(results.input_value, results.input2_value)
	
	# Plot the coverage for busco and cegma id

	if results.step_value == "busco_cegma_cov" and results.input_value and results.input2_value and results.output_value:
		busco_cegma_cov(results.input_value, results.input2_value, results.output_value)
	
	if results.step_value == "random_coverage" and results.input_value and results.output_value:
		random_coverage(results.input_value, results.output_value)
	
	# List scaffolds that have Busco and Cegma ids:
	
	if results.step_value == "scaffolds_busco_cegma" and results.input_value and results.input2_value and results.input3_value and results.output_value:
		scaffolds_busco_cegma(results.input_value, results.input2_value, results.input3_value, results.output_value)

	# Merge reseq data
	
	if results.step_value == "cat_reseq" and results.input_value and results.email_value and results.output_value:
		cat_reseq(results.input_value, results.email_value, results.output_value)
		
	# Map reseq data
	
	if results.step_value == "bwa_reseq" and results.input_value and results.input2_value and results.input3_value and results.output_value:
		bwa_reseq(results.input_value, results.input2_value, results.input3_value, results.output_value)

	# coverage manual and blast manual

	if results.step_value == "coverage_manual" and results.input_value and results.input2_value and results.output_value:
		coverage_manual(results.input_value, results.input2_value, results.output_value)


	if results.step_value == "blast_manual" and results.input_value and results.input2_value and results.output_value:
		blast_manual(results.input_value, results.input2_value, results.output_value)
			
	# size read from fastq

	if results.step_value == "read_size" and results.input_value and results.output_value:
		read_size(results.input_value, results.output_value)					

																									
if __name__ == "__main__":
	main(sys.argv[1:])
