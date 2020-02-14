# 1. data for triangles
# sets of non-rare all homoz. and heteroz.

from collections import defaultdict

for sp in ['1_Tdi', '1_Tps', '2_Tcm', '2_Tsi', '3_Tce', '3_Tms', '4_Tbi', '4_Tte', '5_Tge', '5_Tpa']:
    triangle_data = defaultdict(int)
    variant_file = 'data/' + sp + '.SNP_filter_passed.tsv'
    triangle_file = 'data/' + sp + '_trinalge_SNP_filter_passed.tsv'
    homozygous_variants_file = 'data/' + sp + '_homozygous_SNP_filter_passed.tsv'
    heterozygous_variants_file = 'data/' + sp + '_heterozygous_SNP_filter_passed.tsv'
    with open(variant_file) as f, open(homozygous_variants_file, 'w') as homo_var, open(heterozygous_variants_file, 'w') as het_var:
        for line in f:
            variant = line.split()
            reference = variant[0][:10] + '8' + variant[0][11:]
            heterozygots = sum([v == '0/1' for v in variant])
            homozygots = sum([v == '1/1' for v in variant])
            # save for triangle data to a dict of heterozygots and alleles
            triangle_data[(heterozygots, heterozygots + (2 * homozygots))] += 1
            if heterozygots == 0 and homozygots > 1:
                homo_var.write(reference + '\t' + variant[1] + '\n')
            if homozygots == 0 and heterozygots > 1:
                het_var.write(reference + '\t' + variant[1] + '\n')
    with open(triangle_file, 'w') as triangle_f:
        for het_alleles in triangle_data:
            triangle_f.write('%d\t%d\t%d\n' % (het_alleles[0], het_alleles[1], triangle_data[het_alleles]))
