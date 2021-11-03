import sys
import os
import getopt
import decimal
from decimal import *
import re
import collections

try:
	opts, args = getopt.getopt(sys.argv[1:], 'a:g:o:Dh')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

align_file_name = None
outprefix    = "testout"
genome_prefix = None
no_div_out = False

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h'):
		print("\n**** TE_align_to_gff.py | Written by DJP, 19/10/21 in Python 3.5 in Lausanne, Swiss ****\n")
		print("This code takes Repeat Landscape alignment (algn) files from repeatmasker, and produces a gff file for all alignments with a CpG Kimura distance estimate.")
		print("Note, no strand infomation is outputted in the gff, and the CpG Kimura distances are given in the 6th column.")	
		print("\n**** Usage****\n")
		print("python3 TE_align_to_gff.py -a [alignment file] -g [genome scaffold prefix] -o [output file prefix]\n")
		print("\n**** Options ****\n")
		print("-a\tAlignment file name\n-g\tGenome scaffold prefix - whatever characters start each genome scaffold\n-o\tPrefix for output files\n-D\t[OPTIONAL] Specify if you don't want CpG Kimura distances outputting in the gff file. Default OFF. \n")
		print("\n**** Example ****\n")
		print("python3 TE_align_to_gff.py -a algn_files/1_Tdi_b3v08.fasta.align -g 1_Tdi_b3v08 -o 1_Tdi_b3v08 \n\n")

		sys.exit(2)
	elif opt in ('-a'):
		align_file_name = arg
	elif opt in ('-g'):
		genome_prefix   = arg
	elif opt in ('-D'):
		no_div_out = True
	elif opt in ('-o'):
		outprefix = arg
	else:
		print("i dont know")
		sys.exit(2)


if align_file_name == None:
	print("Error. No alignment file specified.")
	sys.exit(2)

if genome_prefix == None:
	print("Error. No genome prefix set. This is whatever all scaffolds in the genome file begin with.")
	sys.exit(2)
else:
	print("\nThis should be the genome prefix you are using:")
	print(genome_prefix)

print("\n\nNote only records with a CpG Kimura distance estimate will be outputted.\n\n")

align_dict = {}
N_align = 0
N_align_w_div = 0

if no_div_out == False:
	out_file = open(outprefix + "_align.gff", "w")
else:
	print("-D option selected. CpG Kimura distance estimates will not be outputted to the gff file.")
	out_file = open(outprefix + "__no_div_align.gff", "w")

##### read align file

align_file = open(align_file_name)
for line in align_file:
	line_o = line.rstrip("\n")
	line = line.rstrip("\n").split(" ")
	if len(line) > 7:
		if line[0] != "":
			scaf_name = line[4] + "_" + line[5] + "_" + line[6]
			if scaf_name.startswith(genome_prefix):
				scaf_name_w = line[4]
				start_c = line[5] 
				end_c = line[6]
				N_align = N_align + 1
				TE_class = None
				for el in line:
					if "#" in el:
						TE_class = el
				if TE_class == None:
					print("Error getting TE class in align file. Exit")
					sys.exit()
						

	if "Kimura (with divCpGMod)" in line_o:
		#print(line_o)
		
		if no_div_out == False:
			out_file.write(scaf_name_w + "\t" +  "RepeatMasker"  + "\t" + "similarity" + "\t"  + start_c  + "\t" + end_c + "\t" + line_o.split("=")[1].strip() + "\t.\t.\tTarget " + TE_class + " " + scaf_name_w + "_" + start_c  + "_" + end_c + "\n")
		else:
			out_file.write(scaf_name_w + "\t" +  "RepeatMasker"  + "\t" + "similarity" + "\t"  + start_c  + "\t" + end_c + "\t" + "." + "\t.\t.\tTarget " + TE_class + " " + scaf_name_w + "_" + start_c  + "_" + end_c + "\n")
		
		N_align_w_div = N_align_w_div + 1
align_file.close()

print("\nNumber of alignments in " + align_file_name  + " = " + str(N_align))
print("Number of alignments in " + align_file_name  + " with a CpG Kimura distance estimate = " + str(N_align_w_div))

print("\n\nFinished Bero\n\n")

