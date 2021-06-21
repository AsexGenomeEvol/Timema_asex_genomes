import sys
import os
import getopt
import decimal
from decimal import *
import re
import collections

try:
	opts, args = getopt.getopt(sys.argv[1:], 'g:o:f:s:p:G:h')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

gff_file_name = None
BEB_Pr_thresh = decimal.Decimal(0.95)
outprefix    = "testout"
genome_file_name = None

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h'):
		print("\n**** pos_poly.py | Written by DJP, 21/06/21 in Python 3.5 in Lausanne, Swiss ****\n")
		print("Matches polymorphism data to positive selection data.")	
		print("\n**** Usage****\n")
		print("python3 pos_poly.py -g [gff file] -G [genome fasta file] -f [fasta alingments from selectome] -s [vcf file] -p [BEB positively selected sites from selectome] -o [output prefix]\n\n")
		sys.exit(2)
	elif opt in ('-g'):
		gff_file_name  = arg
	elif opt in ('-G'):
		genome_file_name  = arg
	elif opt in ('-f'):
		fasta_align_dir  = arg
	elif opt in ('-s'):
		SNP_file_name  = arg
	elif opt in ('-p'):
		positive_sel_sites_file_name  = arg
	elif opt in ('-o'):
		outprefix = arg
	else:
		print("i dont know")
		sys.exit(2)


def C_DNA(dna):
	dna = dna.upper()
	replacement1 = dna.replace('A', 't')
	replacement2 = replacement1.replace('T', 'a')
	replacement3 = replacement2.replace('C', 'g')
	replacement4 = replacement3.replace('G', 'c')
	return(replacement4.upper())

codon_to_aa_dict = {}

codon_to_aa_dict["TTT"] = "F"
codon_to_aa_dict["TTC"] = "F"
codon_to_aa_dict["TTA"] = "L"
codon_to_aa_dict["TTG"] = "L"
codon_to_aa_dict["TCT"] = "S"
codon_to_aa_dict["TCC"] = "S"
codon_to_aa_dict["TCA"] = "S"
codon_to_aa_dict["TCG"] = "S"
codon_to_aa_dict["TAT"] = "Y"
codon_to_aa_dict["TAC"] = "Y"
codon_to_aa_dict["TAA"] = "STOP"
codon_to_aa_dict["TAG"] = "STOP"
codon_to_aa_dict["TGT"] = "C"
codon_to_aa_dict["TGC"] = "C"
codon_to_aa_dict["TGA"] = "STOP"
codon_to_aa_dict["TGG"] = "W"
codon_to_aa_dict["CTT"] = "L"
codon_to_aa_dict["CTC"] = "L"
codon_to_aa_dict["CTA"] = "L"
codon_to_aa_dict["CTG"] = "L"
codon_to_aa_dict["CCT"] = "P"
codon_to_aa_dict["CCC"] = "P"
codon_to_aa_dict["CCA"] = "P"
codon_to_aa_dict["CCG"] = "P"
codon_to_aa_dict["CAT"] = "H"
codon_to_aa_dict["CAC"] = "H"
codon_to_aa_dict["CAA"] = "Q"
codon_to_aa_dict["CAG"] = "Q"
codon_to_aa_dict["CGT"] = "R"
codon_to_aa_dict["CGC"] = "R"
codon_to_aa_dict["CGA"] = "R"
codon_to_aa_dict["CGG"] = "R"
codon_to_aa_dict["ATT"] = "I"
codon_to_aa_dict["ATC"] = "I"
codon_to_aa_dict["ATA"] = "I"
codon_to_aa_dict["ATG"] = "M"
codon_to_aa_dict["ACT"] = "T"
codon_to_aa_dict["ACC"] = "T"
codon_to_aa_dict["ACA"] = "T"
codon_to_aa_dict["ACG"] = "T"
codon_to_aa_dict["AAT"] = "N"
codon_to_aa_dict["AAC"] = "N"
codon_to_aa_dict["AAA"] = "K"
codon_to_aa_dict["AAG"] = "K"
codon_to_aa_dict["AGT"] = "S"
codon_to_aa_dict["AGC"] = "S"
codon_to_aa_dict["AGA"] = "R"
codon_to_aa_dict["AGG"] = "R"
codon_to_aa_dict["GTT"] = "V"
codon_to_aa_dict["GTC"] = "V"
codon_to_aa_dict["GTA"] = "V"
codon_to_aa_dict["GTG"] = "V"
codon_to_aa_dict["GCT"] = "A"
codon_to_aa_dict["GCC"] = "A"
codon_to_aa_dict["GCA"] = "A"
codon_to_aa_dict["GCG"] = "A"
codon_to_aa_dict["GAT"] = "D"
codon_to_aa_dict["GAC"] = "D"
codon_to_aa_dict["GAA"] = "E"
codon_to_aa_dict["GAG"] = "E"
codon_to_aa_dict["GGT"] = "G"
codon_to_aa_dict["GGC"] = "G"
codon_to_aa_dict["GGA"] = "G"
codon_to_aa_dict["GGG"] = "G"

CDS_scaf_dict   = {}
CDS_coord_dict  = {}
CDS_strand_dict = {}
seen_CDS = set()

species = gff_file_name.strip().split(".gff")[0].split("/")[-1].split("_")[0]
print(species)
### read in file
gff_file = open(gff_file_name)
for line in gff_file:
	line = line.rstrip("\n").split("\t")
	feat = line[2]
	scaf = line[0]
	start_c = line[3]
	end_c   = line[4]
	strand  = line[6] 
	if feat == "CDS":
		ID = line[8].split("Parent=")[1].split(";")[0]
		if strand not in set(["+", "-"]):
			print("strand missing. EXIT")
			sys.exit(2)
		#print(line)
		#print(ID)
		if ID not in seen_CDS:
			seen_CDS.add(ID)
			CDS_scaf_dict[ID] = set([scaf])
			CDS_coord_dict[ID] = [(start_c, end_c)]
			CDS_strand_dict[ID] = set([strand])
		else:
			rec = CDS_scaf_dict.get(ID)
			rec.add(scaf)
			CDS_scaf_dict[ID] = rec
			
			rec_c = CDS_coord_dict.get(ID)
			rec_c.append((start_c, end_c))
			CDS_coord_dict[ID] = rec_c
			
			rec_s = CDS_strand_dict.get(ID)
			rec_s.add(strand)
			CDS_strand_dict[ID] = rec_s
			
gff_file.close()



#### get codon coords 
### check CDSs do not span multiple scafs

gene_codon_pos_to_genome_coord_dict = {}

for el in CDS_scaf_dict:
	rec = CDS_scaf_dict.get(el)
	rec_c = CDS_coord_dict.get(el)
	rec_s = CDS_strand_dict.get(el)
	#print(rec_s)
	if len(rec) != 1:
		print("BAD")
		print(el)
		print(rec)
		sys.exit(2)
		
	full_numbers = []
	
	if list(rec_s)[0] == "+":
		for c in rec_c:
			for n in range(int(c[0]), int(c[1]) + 1):
				full_numbers.append(n)
		full_numbers = sorted(full_numbers)	
	
	elif list(rec_s)[0] == "-":
		for c in rec_c:
			n_sect = []
			for n in range(int(c[0]), int(c[1]) + 1):
				full_numbers.append(n)
			full_numbers = sorted(full_numbers, reverse=True)
	#print(full_numbers)
	## check divides by 3
	if not len(full_numbers) % 3 == 0:
		print("CDS not divisable by 3. EXIT")
		sys.exit(2)
	codon_N = 1
	for f in range(1, len(full_numbers), 3):
		curr_coords = [full_numbers[f -1], full_numbers[f], full_numbers[f +1]] 
		#print(curr_coords)
		gene_codon_pos_to_genome_coord_dict[el + "_codon_" + str(codon_N)]	= (list(rec)[0], curr_coords)
		codon_N = codon_N  + 1


# print(gene_codon_pos_to_genome_coord_dict)
# 		

###############################################################################################
#### alignment codon to CDS codon
#### positive sel results specifiy which codon in the alignment is selected. need to convert to the CDS one calced above.

## read all fasta alignments

### takes dir of fasta files, unwraps them, and adds seqs to a dict with the addition of the file id
def fasta_dir_to_dict(in_seq_dict):
	seq_dict = {}
	path = in_seq_dict
	for path, subdirs, files in os.walk(path):
		for name in files:
			if name.endswith("ORI"):
				#print (os.path.join(path, name))
				
				## file ID ## might need to edit this
				file_id = name.split("/")[-1].split("_to_")[0]
				#print(file_id)
				in_fasta_file_name = os.path.join(path, name)
				output_fasta_name = in_fasta_file_name + ".TEMP_extract_fasta_file" 
				
				output_file = open(output_fasta_name, "w")
				#print("\nUnwrapping fasta file")
				count = 0
				in_file = open(in_fasta_file_name)
				for line in in_file:
					count = count + 1
					line = line.rstrip("\n")
					if line.startswith(">") and count == 1:
						output_file.write(line + "\n")
					elif line.startswith(">") and count > 1:
						output_file.write("\n" + line + "\n")
					else: 
						output_file.write(line)	
				
				output_file.close()
				
				
				### add seqs to dictionary                  	
				name_list = []
				seq_list = []
				
				done = 0
				seq_file_1 = open(output_fasta_name)
				for line in seq_file_1:
					lineA = line.rstrip("\n")
					if lineA.startswith(">"):
						lineB = lineA.lstrip(">")
						name_list.append(lineB)
					else:
						seq_list.append(lineA)
						done = done + 1
						seq_len = len(lineA)
				
				for element in range(0,len(name_list)):
					name1 = name_list[element]
					seq1 = seq_list[element].replace(" ", "") ## remove gaps if seq comes from gblocks 
					seq_dict[name1 + "__FILEID__" + file_id] = seq1
			
				## tidyup
				seq_file_1.close()
				os.remove(output_fasta_name)
	
				#print("Read " + str(done) + " sequences from " + in_fasta_file_name)
	
	return(seq_dict)


align_seq_dict = fasta_dir_to_dict(fasta_align_dir )
print(len(align_seq_dict))

spHOG_to_genename_dict = {}
align_codon_to_CDS_codon_dict = {}
CDS_codon_seq_dict = {}

for gene in align_seq_dict:
	gene_name = gene.split("_HOG_")[0][4:]
	HOG_name  = "HOG_" + gene.split("_HOG_")[1].split("__")[0]
	sp        = gene.split("_")[0]
	
	spHOG_to_genename_dict[sp + "__" + HOG_name] = gene_name
	
	seq = align_seq_dict.get(gene)
	align_codon_N = 0
	CDS_codon_N   = 0
	for i in range(1,len(seq), 3):
		curr_codon = seq[i-1] + seq[i] + seq[i+1]
		align_codon_N = align_codon_N + 1
		if curr_codon != "---":
			CDS_codon_N = CDS_codon_N  + 1
			align_codon_to_CDS_codon_dict[gene_name + "_codon_" + str(align_codon_N)] = CDS_codon_N
			CDS_codon_seq_dict[gene_name + "_codon_" + str(CDS_codon_N)] = curr_codon

		#print(curr_codon + " " + str(align_codon_N) + " " + str(CDS_codon_N))
	
	# print(gene)
	# print(gene_name )
	# print(HOG_name )
	# print(sp)
	#print(seq)

#print(CDS_codon_seq_dict)
#print(align_codon_to_CDS_codon_dict)
#


#### get genome seqs

## takes fasta file, unwraps it, and adds seqs to a dict
def fasta_to_dict(in_fasta_file_name):
	output_fasta_name = in_fasta_file_name + ".TEMP_extract_fasta_file" 
	
	output_file = open(output_fasta_name, "w")
	print("\nUnwrapping fasta file")
	count = 0
	in_file = open(in_fasta_file_name)
	for line in in_file:
		count = count + 1
		line = line.rstrip("\n")
		if line.startswith(">") and count == 1:
			output_file.write(line + "\n")
		elif line.startswith(">") and count > 1:
			output_file.write("\n" + line + "\n")
		else: 
			output_file.write(line)	
	
	output_file.close()
	
	
	### add seqs to dictionary
	name_list = []
	seq_list = []
	seq_dict = {}
	
	done = 0
	seq_file_1 = open(output_fasta_name)
	for line in seq_file_1:
		lineA = line.rstrip("\n")
		if lineA.startswith(">"):
			lineB = lineA.lstrip(">")
			name_list.append(lineB)
		else:
			seq_list.append(lineA)
			done = done + 1
			seq_len = len(lineA)
	
	for element in range(0,len(name_list)):
		name1 = name_list[element]
		seq1 = seq_list[element].replace(" ", "") ## remove gaps if seq comes from gblocks 
		seq_dict[name1] = seq1

	## tidyup
	seq_file_1.close()
	os.remove(output_fasta_name)
	
	print("Read " + str(done) + " sequences from " + in_fasta_file_name)
	
	return(seq_dict)


genome_seq_dict = fasta_to_dict(genome_file_name)


###################################################################################################
### polymorphism data

all_poly_SNPs = set()
poly_dict = {}
SNP_file = open(SNP_file_name)
for line in SNP_file:
	if not line.startswith("#"):
		if "PASS" in line: ## get snps that pass the thresh
			line = line.strip().split("\t")
			scaf = line[0].replace("_b3v06_", "_b3v08_")
			pos    = line[1]
			all_poly_SNPs.add(scaf + "_pos_" + pos)
			poly_dict[scaf + "_pos_" + pos]  = [line[3], line[4]]

# print(poly_dict)
print(len(all_poly_SNPs))


###################################################################################################
### workout if selected codons are polymorphic

positive_sel_sites_file = open(positive_sel_sites_file_name)
trans_sites  = 0
genome_sites = 0
outfile = open(outprefix + "_poly_or_fixed_.csv", "w")
outfile.write("HOG" + "," + "align_pos" +  "," + "gene_name" + "," + "gene_codon" + "," + "genome_coords" + "," + "poly_nucl" + "," + "poly_aa" + "\n")
for line in positive_sel_sites_file:
	line = line.strip().split(",")
	sp = line[3]
	HOG = line[7]
	codon_pos = line[4]
	if sp == species:
		# print(line)
		# print(sp + "\t" + HOG + "\t" + codon_pos)
		gene_name = spHOG_to_genename_dict.get(HOG)
		
		## take only genome genes
		if "TRINITY" in gene_name:
			trans_sites = trans_sites + 1
		else:
			genome_sites = genome_sites + 1
			#print(sp + "\t" + HOG + "\t" + codon_pos)
			# print(gene_name)
			# print(codon_pos)			
			gene_codon    = align_codon_to_CDS_codon_dict.get(gene_name + "_codon_" + codon_pos)
			CDS_codon_seq = CDS_codon_seq_dict.get(gene_name + "_codon_" + str(gene_codon))
			
			genome_coords = gene_codon_pos_to_genome_coord_dict.get(gene_name + "_codon_" + str(gene_codon))
			
			# print(gene_codon)
			#print(genome_coords)
			if genome_coords == None:
				print("Something wrong with the following record. Skipping.")
				print(sp + "\t" + HOG + "\t" + codon_pos)
			else:
				all_genome_coords = [genome_coords[0] +  "_pos_" + str(genome_coords[1][0]), genome_coords[0] +  "_pos_" + str(genome_coords[1][1]), genome_coords[0] +  "_pos_" + str(genome_coords[1][2])]
				genome_seq       = genome_seq_dict.get(genome_coords[0])
				gene_strand      = list(CDS_strand_dict.get(gene_name))[0]
				#print(gene_strand)
				
				if gene_strand == "+":
					genome_seq_codon = genome_seq[genome_coords[1][0] -1]  + genome_seq[genome_coords[1][1] -1] + genome_seq[genome_coords[1][2] -1]
					genome_seq_codon = genome_seq_codon.upper()
				elif gene_strand == "-":
					genome_seq_codon = C_DNA(genome_seq[genome_coords[1][0] -1]  + genome_seq[genome_coords[1][1] -1] + genome_seq[genome_coords[1][2] -1])
				else:
					print("error")
					sys.exit(2)
					
				### check i got the right positions! if not exit.
				if genome_seq_codon != CDS_codon_seq:
					print("CDS codon from alignment and genome codon do not match. Exiting.")
					print(gene_name + "_codon_" + str(gene_codon))
					print(CDS_codon_seq)
					print(all_genome_coords)
					print(genome_seq_codon)
					sys.exit(2)
					
					

				genome_coord_out = ""
				
				#print(all_genome_coords)
				#print(genome_seq_codon)
				
				### make alt codon
				poly_aa   = "fixed"
				poly_nucl = "fixed"
				
				alt_seq = ""
				curr_base = 0
				for c in all_genome_coords:
					genome_coord_out = genome_coord_out + ";" + c
					
					if c in all_poly_SNPs:
						poly_nucl = "poly"
						poly_bases = poly_dict.get(c)
						curr_genome_base = genome_seq[genome_coords[1][curr_base] -1]
						for b in poly_bases:
							if b != curr_genome_base:
								alt_seq = alt_seq + b
								break
						
						#print(poly_bases )
					
					else:
						alt_seq = alt_seq + genome_seq[genome_coords[1][curr_base] -1]

					curr_base = curr_base + 1	

						
				#print(alt_seq)
				
				alt_codon = ""
				if gene_strand == "-":
					alt_codon = C_DNA(alt_seq)
				else:
					alt_codon = alt_seq.upper()
					
				#print(alt_codon)
				
				genome_seq_aa = codon_to_aa_dict.get(genome_seq_codon)
				alt_seq_aa    = codon_to_aa_dict.get(alt_codon)
				
				#print(genome_seq_aa)
				#print(alt_seq_aa)				
				
				if genome_seq_aa != alt_seq_aa:
					poly_aa = "poly" 
				
				genome_coord_out = genome_coord_out.strip(";")
				#print(poly_aa)
				#outfile.write(gene_name + "," + str(gene_codon) + "," + genome_coords[0] +  "_pos_" + str(genome_coords[1][0]) + ";" + genome_coords[0] +  "_pos_" + str(genome_coords[1][1]), + ";" + genome_coords[0] +  "_pos_" + str(genome_coords[1][2]) + "," + poly + "\n")
				outfile.write(HOG + "," + str(codon_pos) +  "," + gene_name + "," + str(gene_codon) + "," + genome_coord_out + "," + poly_nucl + "," + poly_aa + "\n")
				
outfile.close()
print(trans_sites)
print(genome_sites)

print("\n\n\nDone Robin Wednesbury\n\n\n")


sys.exit(1)

