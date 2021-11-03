import sys
import os
import getopt
import decimal
from decimal import *
import re
import collections

try:
	opts, args = getopt.getopt(sys.argv[1:], 'g:c:m:o:h')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

gff_file_name = None
count_file_name = None
min_feat_len = 80
outprefix    = "testout"

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h'):
		print("\n**** TE_coverage_from_align.py | Written by DJP, 20/10/21 in Python 3.5 in Lausanne, Swiss ****\n")
		print("Adds CpG Kimura distance estimates and feature lengths (from the gff produced by TE_align_to_gff.py) to read counts from HTseq.")	
		print("\n**** Usage****\n")
		print("python3 TE_coverage_from_align.py -g [gff file name] -c [HTseq count file name] -o [output file prefix]\n")
		print("\n**** Options ****\n")
		print("-g\tgff file name (the gff produced by TE_align_to_gff.py)\n-c\tHTseq count file name\n-o\tPrefix for output files\n-m\t[OPTIONAL] Minumum feature length. Features shorter than this will be discarded in the output. Default = 80. \n")
		print("\n**** Example ****\n")
		print("python3 TE_coverage_from_align.py -g 1_Tdi_b3v08_align.gff -c 1_Tdi_htseq_unq_withUNIQID_from_align.counts -o 1_Tdi \n\n")		
		sys.exit(2)
	elif opt in ('-g'):
		gff_file_name  = arg
	elif opt in ('-c'):
		count_file_name  = arg
	elif opt in ('-m'):
		min_feat_len  = arg
	elif opt in ('-o'):
		outprefix = arg
	else:
		print("i dont know")
		sys.exit(2)

if gff_file_name == None:
	print("\n\nERROR. No gff file specified\n\n")
	sys.exit(2)

if count_file_name == None:
	print("\n\nERROR. No count file specified\n\n")
	sys.exit(2)

try:
	min_feat_len = int(min_feat_len)
	print("\nMinumum feature length = " + str(min_feat_len) + ". Features shorter than this will be discarded in the output. Change this value with -m ")
except:
	print("\nError: Specified minumum feature length is not an integer. Please specify an integer value.")
	sys.exit(2)
	
### get scaf positon, feat len, and Kimura (with divCpGMod) from gff file (from TE_align_to_gff.py)

gff_dict = {}

gff_file = open(gff_file_name)
for line in gff_file:
	line = line.rstrip("\n").split("\t")
	start_c = line[3]
	end_c   = line[4]
	len_feat = int(end_c) - int(start_c) + 1
	scaf_ID  = line[0] + "_" + start_c + "_" + end_c
	feature_ID = line[8].replace("Target ", "")
	div_k = line[5]
	gff_dict[feature_ID] = [scaf_ID, len_feat, div_k]
	
	# print(line)
	# print(len_feat )
	# print(feature_ID)
	# print(scaf_ID)
	# print(div_k)
	# print("")
gff_file.close()


###### add info to count file

out_file = open(outprefix + "_TE_cov_align.txt", "w")
out_file.write("scaf_ID" + "\t" + "feat_len" + "\t" + "Kimura_with_divCpGMod" + "\t" +  "Nreadcounts" + "\t" + "TE_A_class\tTE_A_class_s" + "\t" + "feat_name" + "\n")

count_file = open(count_file_name)
for line in count_file:
	line = line.rstrip("\n")
	if not line.startswith("__"):
		line = line.split("\t")
		feat = line[0]
		count = line[1]
		gff_rec = gff_dict.get(feat)
		scaf_ID = gff_rec[0]
		feat_len = gff_rec[1]
		div = gff_rec[2]
		# print(scaf_ID)
		# print(feat_len)
		# print(div)
		TE_A_class = feat.split(" ")[0]
		TE_A_class_s = TE_A_class.split("#")[1]
		if int(feat_len) >= min_feat_len:
			out_file.write(scaf_ID + "\t" + str(feat_len) + "\t" + str(div) + "\t" +  count + "\t" + TE_A_class + "\t" + TE_A_class_s + "\t" + feat + "\n")			
	
count_file.close()

print("\n\nFinished Sanda\n\n")











