#! /usr/bin/env python

import os
import sys
from Bio import SeqIO
from optparse import OptionParser

parser = OptionParser(usage='''%prog [-i] [-o] [-s] [-l]
*************************************************************************
This programm sorts a fasta file acording to header (case sensitive)

BioPython needed.
http://biopython.org/wiki/Main_Page
contact: jbast@gwdg.de
*************************************************************************
''', version="%prog 1.0")

parser.add_option("-i", "--infile", dest="InputFile", type='string', help='specify input file', default=None)
parser.add_option("-o", "--outfile", dest="OutputFile", type='string', help='specify output file',default=None)


(options, args) = parser.parse_args()

if options.OutputFile is None:
    parser.error("no output file given")
if options.InputFile is None:
    parser.error("no input file given")

if os.path.exists(options.OutputFile):
	sys.stderr.write('\nFile exists!\n\n')
	overwrite=raw_input('\noverwrite? (y/n): ')
	if overwrite != 'y':
		sys.exit()

#handle = open("seqs.fa", "rU")
with open (options.OutputFile, 'w') as handle:
	#for seq_record in SeqIO.parse(options.InputFile, 'fasta'):
	l = SeqIO.parse(options.InputFile, "fasta")
	sortedList = [f for f in sorted(l, key=lambda x : x.id)]
	for s in sortedList:
   		#print s.description
   		#print str(s.seq)
   		handle.write('>%s\n%s\n' % (s.description, str(s.seq)))

sys.stderr.write('\n written to %s\n\n' % (options.OutputFile))
