#! /usr/bin/env python

import os
import sys
from Bio import SeqIO
from optparse import OptionParser


parser = OptionParser(usage='''%prog [-i] [-o] [-s] [-l]
*************************************************************************
This programm will convert fasta interleaved to fasta non-interleaved
and adds the sequence length
size-threshold optional (-s)
append sequence length (-l)

BioPython needed.
http://biopython.org/wiki/Main_Page
contact: jbast@gwdg.de
*************************************************************************
''', version="%prog 1.0")

parser.add_option("-i", "--infile", dest="InputFile", type='string', help='specify input file', default=None)
parser.add_option("-o", "--outfile", dest="OutputFile", type='string', help='specify output file',default=None)
parser.add_option("-s", "--size", dest="sizeInput", type='int', help='specify minimum contig size',default=0)
parser.add_option("-l", "--length", dest="length", action="store_true", help='append length', default=False)
parser.add_option("-t", "--type", dest="type", type='string', help='specify fasta or fastq',default=None)

(options, args) = parser.parse_args()

if options.OutputFile is None:
    parser.error("no output file given")
if options.InputFile is None:
    parser.error("no input file given")
if options.type is None:
    parser.error("no type given. specify fasta or fastq")

if os.path.exists(options.OutputFile):
	sys.stderr.write('\nFile exists!\n\n')
	overwrite=raw_input('\noverwrite? (y/n): ')
	if overwrite != 'y':
		sys.exit()
		

sys.stderr.write('\nGoing through file...\n\n')	
with open (options.OutputFile, 'w') as handle:

#here, a index will work in ordered  (see cookbook)


#with open (options.OutputFile, 'w') as handle:

	for seq_record in SeqIO.parse(options.InputFile, options.type):
		
		if len(seq_record) >= options.sizeInput:
			if options.length == True:
				handle.write('>%s\t%s\n%s\n' % (seq_record.id, len(seq_record), seq_record.seq))
			else:
				handle.write('>%s\n%s\n' % (seq_record.id, seq_record.seq))

		
sys.stderr.write('\n written to %s\n\n' % (options.OutputFile))
