#!/usr/bin/env python3

from collections import defaultdict
import argparse
from sys import stdout
from sys import stderr
import numpy as np
import gzip

class Genotype:
    def __init__(self, genotype_string):
        parsed_genotype = genotype_string.split(':')
        self.gt = parsed_genotype[0]
        self.gq = parsed_genotype[3]

    def __str__(self):
        return(self.gt + '\t' + self.gq)

    def __repr__(self):
        return(self.gt + '\t' + self.gq)

class Variant:
    def __init__(self, line, ind, ind_order):
        parsed_line = line.rstrip('\n').split('\t')
        sorted_indexes_of_individuals = [9 + i[0] for i in sorted(enumerate(ind_order), key=lambda x:x[1])]
        self.line = parsed_line
        self.scf = parsed_line[0]
        self.pos = int(parsed_line[1])
        info_dict = dict(s.split('=', 1) for s in parsed_line[7].split(';') if "=" in s)
        self.type = info_dict['SVTYPE']
        self.len = int(info_dict['SVLEN'])
        self.orig = ind
        self.called = parsed_line[9].split(':')[1]
        self.genotypes = [Genotype(parsed_line[i]) for i in sorted_indexes_of_individuals]

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

    def vcfPrint(self):
        return("\t".join(self.line) + "\n")

    def isValid(self):
        return(any([gt.gq == "PASS" for gt in self.genotypes]))

#############
# FUNCTIONS #
#############
# I guess I could have made a class for a buffer of overlapping variants
# but I have decided to do it using a few functions

def processes_buffer(overlap_buffer):
    # decide which in overlap_buffer to be merged and if to print them to
    while overlap_buffer:
        to_merge = [0]
        for i in range(1, len(overlap_buffer)):
            # merging conditions
            if variant_agreement(overlap_buffer[0], overlap_buffer[1]):
                to_merge.append(i)

        # positions = np.mean([overlap_buffer[i].pos for i in to_merge])
        final_variant = merge_variants([overlap_buffer[i] for i in to_merge])
        if final_variant.isValid():
            stdout.write(final_variant.vcfPrint())
        else:
            stderr.write(final_variant.vcfPrint())

        # remove merged variants from the list of variants
        for to_del in sorted(to_merge, reverse=True):
            del overlap_buffer[to_del]

def merge_variants(var_list):
    # TODO: merge them in an inteligent manner
    return(var_list[0])

def variant_agreement(var1, var2):
    # TODO: identify if the variants have compatible calls
    # for i in range(6):
    #   comare individual i in var1 and var2
    return(True)

def variant_overlap(var1, var2):
    return(var1.type == var2.type and var1.pos < (var2.pos + abs(var2.len)) and var2.pos < (var1.pos + abs(var1.len)))

##########
# SCRIPT #
##########

parser = argparse.ArgumentParser(description='Reciprocal BLAST: from blast std output 6 generates two tables of orthologs as two with reciprocal hits fulfilling specified criteria')
parser.add_argument('sp', help='The species code (e.g. 3_Tms)')

args = parser.parse_args()
# args.sp = '3_Tms'

stderr.write('# processing ' + args.sp + '\n')
stderr.write('# output will be streamed to stdout\n')
stderr.write('# log of the computation is streamed to stderr and starting with #\n')
stderr.write('# the non # lines on the stderr are filtered variants\n')
stderr.write('# python3 E_structural_variation_calling/merging_paragraph_calls.py 3_Tms > merged_variants.vcf 2> log_and_filtered_variants.vcf\n')

# files = ['data/genotyping/' + args.sp + '_ind0' + str(i) + '_genotyping/genotypes.vcf.gz' for i in range(6)]
# ind = 2
# file = files[ind]
scf2variants = defaultdict(list)

for ind in range(6):
    file = 'data/genotyping/' + args.sp + '_ind0' + str(ind) + '_genotyping/genotypes.vcf.gz'
    stderr.write('# reading ' + file + '\n')
    with gzip.open(file, mode='rt') as vcf_file:
        for line in vcf_file:
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                header = line.rstrip('\n').split('\t')
                ind_order = header[9:15]
                stderr.write('# individual order ' + '\t'.join(ind_order) + '\n')
                continue
            processed_variant = Variant(line, ind, ind_order)
            scf2variants[processed_variant.scf].append(processed_variant)

stderr.write('# variants loaded...\n')

for scf in scf2variants.keys():
# scf = '3_Tms_b3v08_scaf000008'
    scf2variants[scf].sort()
    overlap_buffer = []
    buffer = 1
    for var in scf2variants[scf]:
        if not overlap_buffer:
            overlap_buffer.append(var)
            continue
        if any([ variant_overlap(var, var2) for var2 in overlap_buffer ]):
            overlap_buffer.append(var)
        else:
            if len(overlap_buffer) > 6:
                stderr.write('#! found a buffer of size:' + str(len(overlap_buffer)) + '\n')
                stderr.write('#! ' + str(overlap_buffer[0]) + '\n')
                #break
            if buffer == 3:
                break

            processes_buffer(overlap_buffer)
            buffer += 1
            overlap_buffer = []
            overlap_buffer.append(var)

stderr.write('# Done\n')

# merging only overlaps with complete agreement

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