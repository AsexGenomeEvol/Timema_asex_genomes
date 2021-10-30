# -*- coding: utf-8 -*-

import os
import argparse
import sys
import glob
import os.path
import shutil
import re
from collections import OrderedDict
from difflib import SequenceMatcher
import ast

#### HOW TO USE

## python $SC/te_landscape.py -s repeatmasker2r -i1 1_Tdi_b3v08.fa_mod.html -o 1_Tdi_b3v08_R_data.txt
	
def repeatmasker2r(repeatmasker_file, output_file):

	"""	
	Convert Repeatmasker output to table suitable for R
	"""
		
	if os.path.isfile(output_file):
		os.remove(output_file)
		
	
	script = open(output_file,"a")		
	script.write("Type\tKimura\tGenome_percent\n")				
	script.close()
	
	
	type_TE = []
	line_no_space_list = []		
							
	## Parsing the RepeatMasker file
	
	repeatmasker = open(repeatmasker_file,"r")
	
	canPrintLines = False 
	
	for line in repeatmasker:
		
		if re.search("data.addColumn", line) and re.search("number", line):
			type_TE.append(re.split("[\r\t\n']",line)[3])

		# Insert in a list all lines between data.addRows and varpieData
		
		
		line_no_space = line.replace(" ","")
		if line_no_space.startswith("data.addRows(["):
			canPrintLines = True 
		elif line_no_space.startswith("varpieData"):
			canPrintLines = False 

		if canPrintLines:
			line_no_space_list.append(line_no_space.split(",\n"))		
	#print type_TE
	
	repeatmasker.close()

	# https://stackoverflow.com/questions/1894269/convert-string-representation-of-list-to-list
	kimura_list = []
	
	kimura_index = 1
	while kimura_index < len(line_no_space_list)-1:
		kimura_data_to_list = ast.literal_eval(line_no_space_list[kimura_index][0])
		#print kimura_data_to_list
		
		kimura_data_list_formated = []
		
		for percent in kimura_data_to_list:
			kimura_data_list_formated.append(str(percent))
			 
		kimura_list.append(kimura_data_list_formated)
		#v = [n.strip() for n in x]
		#print v, len(v)
		kimura_index += 1

	#print kimura_list
	
	## Reverse the list
	
	# Before: ['Other', 'DNA/Academ', ... ,'LINE/R2', 'SINE']
	# After: ['SINE', 'LINE/R2', ... , 'DNA/Academ', 'Other']
	
	type_TE = type_TE[::-1]
	
	#print type_TE, len(type_TE)
	
	## Same for kimura list
	
	kimura_list_reverse = []
	
	for kimura in kimura_list:
		#print kimura, len(kimura)
		
		# Reverse all except the 1st value (kimura substitution level)
		kimura_list_reverse.append([kimura[0]] + kimura[1:][::-1])
	
	for kimura in kimura_list_reverse:
	
		#print kimura_list_reverse
		dic_kimura = {}
		data = 1
	
		while data < len(kimura):
		
			dic_kimura[type_TE[data-1]] = kimura[data]
			
			data += 1
	
		#print kimura
		# ['50', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0']
		
		# -> 1 st column is substitution level
		
		script = open(output_file,"a")
		
		if "SINE" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("SINE", kimura[0], dic_kimura["SINE"]))
		else:
			script.write("{}\t{}\t{}\n".format("SINE", kimura[0], "0"))

		if "LINE/Crack" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("LINE/Crack", kimura[0], dic_kimura["LINE/Crack"]))
		else:
			script.write("{}\t{}\t{}\n".format("LINE/Crack", kimura[0], "0"))

		if "LINE/Tx1" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("LINE/Tx1", kimura[0], dic_kimura["LINE/Tx1"]))
		else:
			script.write("{}\t{}\t{}\n".format("LINE/Tx1", kimura[0], "0"))

		if "LINE/Poseidon" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("LINE/Poseidon", kimura[0], dic_kimura["LINE/Poseidon"]))
		else:
			script.write("{}\t{}\t{}\n".format("LINE/Poseidon", kimura[0], "0"))

		if "LINE/R2" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("LINE/R2", kimura[0], dic_kimura["LINE/R2"]))
		else:
			script.write("{}\t{}\t{}\n".format("LINE/R2", kimura[0], "0"))

		if "LINE/Jockey-I" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("LINE/Jockey-I", kimura[0], dic_kimura["LINE/Jockey-I"]))
		else:
			script.write("{}\t{}\t{}\n".format("LINE/Jockey-I", kimura[0], "0"))

		if "LINE/R1" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("LINE/R1", kimura[0], dic_kimura["LINE/R1"]))
		else:
			script.write("{}\t{}\t{}\n".format("LINE/R1", kimura[0], "0"))

		if "LINE/LOA" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("LINE/LOA", kimura[0], dic_kimura["LINE/LOA"]))
		else:
			script.write("{}\t{}\t{}\n".format("LINE/LOA", kimura[0], "0"))

		if "LINE/L2" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("LINE/L2", kimura[0], dic_kimura["LINE/L2"]))
		else:
			script.write("{}\t{}\t{}\n".format("LINE/L2", kimura[0], "0"))

		if "LINE/CR1" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("LINE/CR1", kimura[0], dic_kimura["LINE/CR1"]))
		else:
			script.write("{}\t{}\t{}\n".format("LINE/CR1", kimura[0], "0"))

		if "LINE/RTE" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("LINE/RTE", kimura[0], dic_kimura["LINE/RTE"]))
		else:
			script.write("{}\t{}\t{}\n".format("LINE/RTE", kimura[0], "0"))

		if "LINE" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("LINE", kimura[0], dic_kimura["LINE"]))
		else:
			script.write("{}\t{}\t{}\n".format("LINE", kimura[0], "0"))

		if "LINE/L1" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("LINE/L1", kimura[0], dic_kimura["LINE/L1"]))
		else:
			script.write("{}\t{}\t{}\n".format("LINE/L1", kimura[0], "0"))

		if "LTR/Penelope" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("LTR/Penelope", kimura[0], dic_kimura["LTR/Penelope"]))
		else:
			script.write("{}\t{}\t{}\n".format("LTR/Penelope", kimura[0], "0"))

		if "LTR" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("LTR", kimura[0], dic_kimura["LTR"]))
		else:
			script.write("{}\t{}\t{}\n".format("LTR", kimura[0], "0"))

		if "LTR/Gypsy" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("LTR/Gypsy", kimura[0], dic_kimura["LTR/Gypsy"]))
		else:
			script.write("{}\t{}\t{}\n".format("LTR/Gypsy", kimura[0], "0"))

		if "LTR/Copia" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("LTR/Copia", kimura[0], dic_kimura["LTR/Copia"]))
		else:
			script.write("{}\t{}\t{}\n".format("LTR/Copia", kimura[0], "0"))

		if "LTR/BEL" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("LTR/BEL", kimura[0], dic_kimura["LTR/BEL"]))
		else:
			script.write("{}\t{}\t{}\n".format("LTR/BEL", kimura[0], "0"))

		if "RC/Helitron" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("RC/Helitron", kimura[0], dic_kimura["RC/Helitron"]))
		else:
			script.write("{}\t{}\t{}\n".format("RC/Helitron", kimura[0], "0"))

		if "DNA/ISL2EU" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("DNA/ISL2EU", kimura[0], dic_kimura["DNA/ISL2EU"]))
		else:
			script.write("{}\t{}\t{}\n".format("DNA/ISL2EU", kimura[0], "0"))

		if "DNA/MITE" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("DNA/MITE", kimura[0], dic_kimura["DNA/MITE"]))
		else:
			script.write("{}\t{}\t{}\n".format("DNA/MITE", kimura[0], "0"))

		if "DNA/MuDr" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("DNA/MuDr", kimura[0], dic_kimura["DNA/MuDr"]))
		else:
			script.write("{}\t{}\t{}\n".format("DNA/MuDr", kimura[0], "0"))

		if "DNA/Transib" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("DNA/Transib", kimura[0], dic_kimura["DNA/Transib"]))
		else:
			script.write("{}\t{}\t{}\n".format("DNA/Transib", kimura[0], "0"))

		if "DNA/TcMar" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("DNA/TcMar", kimura[0], dic_kimura["DNA/TcMar"]))
		else:
			script.write("{}\t{}\t{}\n".format("DNA/TcMar", kimura[0], "0"))

		if "DNA/Sola" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("DNA/Sola", kimura[0], dic_kimura["DNA/Sola"]))
		else:
			script.write("{}\t{}\t{}\n".format("DNA/Sola", kimura[0], "0"))

		if "DNA/PiggyBac" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("DNA/PiggyBac", kimura[0], dic_kimura["DNA/PiggyBac"]))
		else:
			script.write("{}\t{}\t{}\n".format("DNA/PiggyBac", kimura[0], "0"))

		if "DNA/P" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("DNA/P", kimura[0], dic_kimura["DNA/P"]))
		else:
			script.write("{}\t{}\t{}\n".format("DNA/P", kimura[0], "0"))

		if "DNA" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("DNA", kimura[0], dic_kimura["DNA"]))
		else:
			script.write("{}\t{}\t{}\n".format("DNA", kimura[0], "0"))

		if "DNA/Maverick" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("DNA/Maverick", kimura[0], dic_kimura["DNA/Maverick"]))
		else:
			script.write("{}\t{}\t{}\n".format("DNA/Maverick", kimura[0], "0"))

		if "DNA/Kolobok" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("DNA/Kolobok", kimura[0], dic_kimura["DNA/Kolobok"]))
		else:
			script.write("{}\t{}\t{}\n".format("DNA/Kolobok", kimura[0], "0"))

		if "DNA/hAT" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("DNA/hAT", kimura[0], dic_kimura["DNA/hAT"]))
		else:
			script.write("{}\t{}\t{}\n".format("DNA/hAT", kimura[0], "0"))

		if "DNA/Harbinger" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("DNA/Harbinger", kimura[0], dic_kimura["DNA/Harbinger"]))
		else:
			script.write("{}\t{}\t{}\n".format("DNA/Harbinger", kimura[0], "0"))

		if "DNA/Ginger" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("DNA/Ginger", kimura[0], dic_kimura["DNA/Ginger"]))
		else:
			script.write("{}\t{}\t{}\n".format("DNA/Ginger", kimura[0], "0"))

		if "DNA/Chapaev" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("DNA/Chapaev", kimura[0], dic_kimura["DNA/Chapaev"]))
		else:
			script.write("{}\t{}\t{}\n".format("DNA/Chapaev", kimura[0], "0"))

		if "DNA/Academ" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("DNA/Academ", kimura[0], dic_kimura["DNA/Academ"]))
		else:
			script.write("{}\t{}\t{}\n".format("DNA/Academ", kimura[0], "0"))

		if "Other" in dic_kimura:
			script.write("{}\t{}\t{}\n".format("Other", kimura[0], dic_kimura["Other"]))
		else:
			script.write("{}\t{}\t{}\n".format("Other", kimura[0], "0"))
	
				
		script.close()
			
	print "Table has been written."

												
def main(argv):
	
	mod=[]
	mod.append('\n%(prog)s -s repeatmasker2r -i1 <hog_file> -i2 <table_file> -o <output_file>')
	
	parser = argparse.ArgumentParser(prog = 'te.py',
                                 usage = "\n".join(mod))

	parser.add_argument('-s', action='store', dest='step_value',
	                    help='Step')
	                                                 	
	parser.add_argument('-i1', action='store', dest='input_value',
	                    help='Input 1')

	parser.add_argument('-i2', action='store', dest='input2_value',
	                    help='Input 2')

	parser.add_argument('-i3', action='store', dest='input3_value',
	                    help='Input 3')
	                    	                    	                    	
	parser.add_argument('-o', action='store', dest='output_value',
	                    help='Output')
	                    	                    	
	parser.add_argument('--version', action='version', version='%(prog)s 1.1')

	results = parser.parse_args()
		
	if results.step_value == "repeatmasker2r" and results.input_value and results.output_value:
		repeatmasker2r(results.input_value, results.output_value)
		
						
if __name__ == "__main__":
	main(sys.argv[1:])
