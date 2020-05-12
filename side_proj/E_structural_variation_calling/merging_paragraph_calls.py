#!/usr/bin/env python3

from collections import defaultdict
import argparse
from sys import stdout
import gzip


parser = argparse.ArgumentParser(description='Reciprocal BLAST: from blast std output 6 generates two tables of orthologs as two with reciprocal hits fulfilling specified criteria')
parser.add_argument('sp', help='The species code (e.g. 3_Tms)')

args = parser.parse_args()
args.sp = '3_Tms'

class genotype:
    def __init__(self, genotype_string):
        parsed_genotype = genotype_string.split(':')
        self.gt = parsed_genotype[0]
        self.filter = parsed_genotype[3]

    def __str__(self):
        return(self.gt + '\t' + self.filter)

    def __repr__(self):
        return(self.gt + '\t' + self.filter)

class variant:
    def __init__(self, line, ind, ind_order):
        parsed_line = line.rstrip('\n').split('\t')
        sorted_indexes_of_individuals = [9 + i[0] for i in sorted(enumerate(ind_order), key=lambda x:x[1])]
        self.line = parsed_line
        self.scf = parsed_line[0]
        self.pos = int(parsed_line[1])
        self.type = parsed_line[5]
        self.orig = ind
        self.called = parsed_line[9].split(':')[1]
        self.genotypes = [genotype(parsed_line[i]) for i in sorted_indexes_of_individuals]

    def __str__(self):
        return(self.scf + "\t" + str(self.pos) + "\t" + self.type)

    def __repr__(self):
        return(self.scf + "\t" + str(self.pos) + "\t" + self.type)

    def __getitem__(self, index):
        return(self.genotypes[index])

    def __lt__(self, other):
         return self.pos < other.pos

    def isInAgreement(self):
        return(self.called == self.genotypes[ind].gt)

    def isPassed(self):
        return(self.called == self.genotypes[ind].gt)

# files = ['data/genotyping/' + args.sp + '_ind0' + str(i) + '_genotyping/genotypes.vcf.gz' for i in range(6)]
# ind = 2
# file = files[ind]
scf2variants = defaultdict(list)

for ind in range(6):
    file = 'data/genotyping/' + args.sp + '_ind0' + str(ind) + '_genotyping/genotypes.vcf.gz'
    with gzip.open(file, mode='rt') as vcf_file:
        for line in vcf_file:
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                header = line.rstrip('\n').split('\t')
                ind_order = header[9:15]
                print(ind_order)
                continue
            processed_variant = variant(line, ind, ind_order)
            scf2variants[processed_variant.scf].append(processed_variant)

# '3_Tms_b3v08_scaf000008' my testing scf

scf = '3_Tms_b3v08_scaf000008'
scf2variants[scf].sort()

overlap_buffer = []

def print_buffer(overlap_buffer):
    # decide which in overlap_buffer to be merged and if to print them to 
    print('Kobyla ma maly bok')
    return(0)

# for scf in scf2variants.keys():
#     scf2variants[scf].sort()
#     overlap_buffer = []
#     for variant in scf2variants[scf]:
# ...
for overlap_buffer in scf2variants[scf]:
    if not overlap_buffer:
        overlap_buffer.append()
        continue
    if any([ variant.overlaps_with(x) for x in overlap_buffer ]):
        overlap_buffer.append()
    else:
        print_buffer(overlap_buffer)
        overlap_buffer = [overlap_buffer]

# fist go, keeping only PASS variants
# merging only overlaps with complete agreement

# GT:
# OLD_GT:
# DP:
# FT:
# AD:ADF:ADR:
# PL:
# GQ:
# PR:
# SR:
# -> seems that nothing is really preserved, the extraction must be done by tag deifnition in every genotype line
# FORMAT/GT: genotype computed by grmpy. Note that the SAMPLE column in the VCF file must match a sample id from the manifest shown above.
# FORMAT/AD, FORMAT/ADF, FORMAT/ADR: depth for each allele (including the reference), total and by strand.
# FORMAT/DP: Total depth used to genotype.
# FORMAT/FT: grmpy filter status for each call, in each sample.
# FORMAT/PL: phred-scaled genotype likelihood.
# GT: 0/0
# OLD_GT: 0/1
# DP: 52
# FT: PASS
# AD: 170,0
# ADF: 101,0
# ADR: 69,0
# PL: 0,117,539
# GQ: 64
# PR: 38,11