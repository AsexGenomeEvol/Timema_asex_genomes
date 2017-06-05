import sys
import os
from Bio import SeqIO

# args :
# SP (1_Tdi)
# REF (b3v04)
# WINDOW (500000)

# prints :
# the ch smaller than half of the window

troot = os.environ['TROOT']
species = sys.argv[1]
version = sys.argv[2]
window_size = int(sys.argv[3])

# $WINDOW < 0.50
treshold = window_size * 0.5
ch = 1

reference_file = troot + '/data/' + species + '/reference/' + species + '_' + version + '.fa'
reference = SeqIO.parse(reference_file, "fasta")

for seq_record in reference:
    if len(seq_record) < treshold:
        sys.stdout.write(ch)
        break
    ch += 1
