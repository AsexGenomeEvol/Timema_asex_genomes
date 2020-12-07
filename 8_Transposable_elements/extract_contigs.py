#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import re
from Bio import SeqIO
from Bio.SeqUtils import GC
from optparse import OptionParser

parser = OptionParser(usage='''%prog [-i] [-o] [-l] [-r] [-h]

*************************************************************************
Extracts contigs from a fasta file using a contig ID list (header part).
Possible to extract all from list or not in list (reverse -r).
-i input (fasta file)
-l list (ID list)
-o output

BioPython needed!
http://biopython.org/wiki/Main_Page
contact: jbast@gwdg.de
*************************************************************************

''', version="%prog 1.0")
parser.add_option("-i", "--input", dest="InputFile", type='string', help='specify file to be searched', default=None)
parser.add_option("-l", "--list", dest="ListFile", type='string', help='search list',default=None)
parser.add_option("-o", "--outfile", dest="OutputFile", type='string', help='specify output file',default=None)
parser.add_option("-r", "--reverse", dest="reverse", action="store_true", help='reverse hit', default=False)

(options, args) = parser.parse_args()
	
if options.OutputFile is None:
	parser.error("no output file given")
if options.InputFile is None:
	parser.error("no input file given")
    
if os.path.exists(options.OutputFile):
	sys.stderr.write('\nFile exists!\n\n')
	overwrite=input('overwrite? (y/n): ')
	if overwrite != 'y':
		sys.exit()

#here, (with open('file', 'w' as handle) you do not need to close the outfile
sys.stderr.write('\nGoing through file...\n\n')


wanted = set()
with open(options.ListFile) as f:
	for line in f:
		line = line.strip()
		if line != "":
           		 wanted.add(line)
fasta_sequences = SeqIO.parse(open(options.InputFile),'fasta')
end = False
with open(options.OutputFile, "w") as f:
	for seq in fasta_sequences:
		if options.reverse == True:
			if seq.id not in wanted:
				SeqIO.write([seq], f, "fasta")
		else:
			if seq.id in wanted:
				SeqIO.write([seq], f, "fasta")

sys.stderr.write('\nDone\n\n')

