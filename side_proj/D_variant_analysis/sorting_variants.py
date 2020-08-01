# 1. data for triangles
# sets of non-rare all homoz. and heteroz.

from collections import defaultdict

for sp in ['1_Tdi', '1_Tps', '2_Tcm', '2_Tsi', '3_Tce', '3_Tms', '4_Tbi', '4_Tte', '5_Tge', '5_Tpa']:
    triangle_data = defaultdict(int)
    variant_file = 'data/SNP_calls/' + sp + '.SNP_filter_passed.tsv'
    triangle_file = 'data/SNP_calls/' + sp + '_trinalge_SNP_filter_passed.tsv'
    homozygous_variants_file = 'data/SNP_calls/' + sp + '_homozygous_SNP_filter_passed.tsv'
    heterozygous_variants_file = 'data/SNP_calls/' + sp + '_heterozygous_SNP_filter_passed.tsv'
    reduced_variants_file = 'data/SNP_calls/' + sp + '_reduced_variants.tsv'
    homo_depth_file = 'data/SNP_calls/' + sp + '_homo_depths.tsv'
    hetero_depth_file = 'data/SNP_calls/' + sp + '_hetero_depths.tsv'
    with open(variant_file) as f, open(homozygous_variants_file, 'w') as homo_var, open(heterozygous_variants_file, 'w') as het_var, open(reduced_variants_file, 'w') as reduced_variants, open(homo_depth_file, 'w') as homo_depth, open(hetero_depth_file, 'w') as hetero_depth:
        for line in f:
            variant = line.split()
            genotypes = [v.split(':')[0] for v in variant[3:]]
            depths = [v.split(':')[2] for v in variant[3:]]
            reference = variant[0][:10] + '8' + variant[0][11:]
            pos = variant[1]
            qual = variant[2]
            heterozygots = sum([g == '0/1' for g in genotypes])
            homozygots = sum([g == '1/1' for g in genotypes])
            # save for triangle data to a dict of heterozygots and alleles
            triangle_data[(heterozygots, heterozygots + (2 * homozygots))] += 1
            if heterozygots == 0 and homozygots > 1:
                homo_var.write(reference + '\t' + pos + '\n')
                homo_depth.write('\t'.join(depths) + '\n')
            if homozygots == 0 and heterozygots > 1:
                het_var.write(reference + '\t' + pos + '\n')
                hetero_depth.write('\t'.join(depths) + '\n')
            reduced_variants.write(reference + '\t' + pos + '\t' + qual + '\t' + '\t'.join(genotypes) + '\t' + '\t'.join(depths) + '\n')
    with open(triangle_file, 'w') as triangle_f:
        for het_alleles in triangle_data:
            triangle_f.write('%d\t%d\t%d\n' % (het_alleles[0], het_alleles[1], triangle_data[het_alleles]))
