import sys
import os
from Bio import SeqIO

# args :
# SP (1_Tdi)
# REF (b3v04)
# WINDOW (500000)

# prints :
# ch from to # (because atlas "to" does not use anyway)
# ch (from + w) (to + w) # (if ch_len - (from + w) > 0.5 * w)
# ch (from + 2 * w) (to + 2 * w) # (if ch_len - (from + 2 * w) > 0.5 * w)
# ...

troot = os.environ['TROOT']
species = sys.argv[1]
version = sys.argv[2]
window_size = int(sys.argv[3])
# $WINDOW < 0.50
treshold = window_size * 0.5

reference_file = troot + '/data/' + species + '/reference/' + species + '_' + version + '.fa'
reference = SeqIO.parse(reference_file, "fasta")

for seq_record in reference:
    scf_len = len(seq_record)
    window_from = 0
    window_to = min(scf_len, window_size - 1)
    while((window_to - window_from) > treshold):
        print(seq_record.name, window_from, window_to, sep = '\t')
        window_from += window_size
        window_to = min(scf_len, window_to + window_size)
    if scf_len < treshold:
        break
