# 1. data for triangles
# sets of non-rare all homoz. and heteroz.

def test_depth(depth):
    if depth == '.':
        return False
    else:
        return int(depth) > 15

from collections import defaultdict
from interval import interval
from sys import stderr

class reference_block:
    def __init__(self, line_tab):
        self.query_from = int(line_tab[3])
        self.query_to = int(line_tab[4])
        self.scf = line_tab[5]
        self.scf_from = int(line_tab[6])
        self.scf_to = int(line_tab[7])
        self.lg = line_tab[8]
        self.lg_from = int(line_tab[9])
        self.lg_to = int(line_tab[10])
    def lg_pos(self, pos):
        aln_size = self.query_to - self.query_from
        ref_aln_size = self.lg_to - self.lg_from
        rel_distance_from_left = (pos - self.query_from) / aln_size
        if self.lg_to > self.lg_from:
            pos = self.lg_from + round(rel_distance_from_left * ref_aln_size)
        else:
            pos = self.lg_to + round((1 - rel_distance_from_left) * ref_aln_size)
        return(pos)
    def scf_pos(self, pos):
        aln_size = self.query_to - self.query_from
        ref_aln_size = self.scf_to - self.scf_from
        rel_distance_from_left = (pos - self.query_from) / aln_size
        if self.scf_to > self.scf_from:
            pos = self.scf_from + round(rel_distance_from_left * ref_aln_size)
        else:
            pos = self.scf_to + round((1 - rel_distance_from_left) * ref_aln_size)
        return(pos)



for sp in ['1_Tdi', '1_Tps', '2_Tcm', '2_Tsi', '3_Tce', '3_Tms', '4_Tbi', '4_Tte', '5_Tge', '5_Tpa']:
    stderr.write("processing {}\n".format(sp))

    variant_file = 'data/SNP_calls/' + sp + '.SNP_filter_passed.tsv'
    mapping_file = 'data/b3v08_anchoring_to_LGs/' + sp + '_scf_block_alignment.tsv'
    reduced_variants_file = 'data/SNP_calls/' + sp + '_reduced_variants.tsv'
    filtered_variants_file = 'data/SNP_calls/' + sp + '_reduced_filtered_variants.tsv'

    # add here also the LG / pos information
    scf2intervals = defaultdict(list)
    scf2mapping = defaultdict(list)
    with open(mapping_file) as mf:
        header = mf.readline()
        for line in mf:
            mapped_block_tab = line.rstrip('\n').split('\t')
            # equvalent of R's mapping_table[abs(mapping_table$block_r_end - mapping_table$block_r_start) / mapping_table$block_size > 0.5, ]
            block_size = int(mapped_block_tab[2])
            block_r_start = int(mapped_block_tab[6])
            block_r_end = int(mapped_block_tab[7])
            if (abs(block_r_end - block_r_start) / block_size) > 0.5:
                mapped_interval = interval([int(mapped_block_tab[3]), int(mapped_block_tab[4])])
                scf2intervals[mapped_block_tab[0]].append(mapped_interval)
                scf2mapping[mapped_block_tab[0]].append(reference_block(mapped_block_tab))

    stderr.write("\tmapping loaded\n")


    with open(variant_file) as vf, open(reduced_variants_file, 'w') as reduced_variants, open(filtered_variants_file, 'w') as filtered_variants:
        unmapped = 0
        filtered = 0
        total = 0
        for line in vf:
            total += 1
            variant = line.split()
            genotypes = [v.split(':')[0] for v in variant[3:]]
            depths = [v.split(':')[2] for v in variant[3:]]
            reference = variant[0][:10] + '8' + variant[0][11:]
            pos = int(variant[1])
            qual = variant[2]
            reduced_variants.write(reference + '\t' + str(pos) + '\t' + qual + '\t' + '\t'.join(genotypes) + '\t' + '\t'.join(depths) + '\n')
            if float(qual) > 300 and any([test_depth(d) for d in depths]):
                mapped_to = list()
                for i, genomic_interval in enumerate(scf2intervals[reference]):
                    if pos in genomic_interval:
                        mapped_to.append(i)
                if len(mapped_to) == 0:
                    # not mapped
                    last_four_cols = "\tNA\tNA\tNA\tNA\n" # these two look the same
                    unmapped += 1

                if len(mapped_to) > 1:
                    # multimapping
                    lg = scf2mapping[reference][mapped_to[0]].lg
                    scf = scf2mapping[reference][mapped_to[0]].scf
                    if not all([scf2mapping[reference][i].lg == lg for i in mapped_to]):
                        lg = "NA"
                    if not all([scf2mapping[reference][i].scf == scf for i in mapped_to]):
                        scf = "NA"
                    last_four_cols = "\t{}\tNA\t{}\tNA\n".format(scf, lg)

                if len(mapped_to) == 1:
                    # unique mapping
                    mapping = scf2mapping[reference][mapped_to[0]]
                    last_four_cols = "\t{}\t{}\t{}\t{}\n".format(mapping.scf, mapping.scf_pos(pos), mapping.lg, mapping.lg_pos(pos))

                filtered_variants.write(reference + '\t' + str(pos) + '\t' + qual + '\t' + '\t'.join(genotypes) + '\t' + '\t'.join(depths) + last_four_cols)

            else:
                filtered += 1

        stderr.write("\tprocessed total {} variants\n".format(total))
        stderr.write("\tfiltered {} variants\n".format(filtered))
        stderr.write("\tunmapped {} variants (not counting multimapping)\n".format(unmapped))
        stderr.write("\tdone\n")
