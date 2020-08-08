# 1. data for triangles
# sets of non-rare all homoz. and heteroz.

def test_depth(depth):
    if depth == '.':
        return False
    else:
        return int(depth) > 15

from collections import defaultdict

for sp in ['1_Tdi', '1_Tps', '2_Tcm', '2_Tsi', '3_Tce', '3_Tms', '4_Tbi', '4_Tte', '5_Tge', '5_Tpa']:
    variant_file = 'data/SNP_calls/' + sp + '.SNP_filter_passed.tsv'
    reduced_variants_file = 'data/SNP_calls/' + sp + '_reduced_variants.tsv'
    filtered_variants_file = 'data/SNP_calls/' + sp + '_reduced_filtered_variants.tsv'
    with open(variant_file) as f, open(reduced_variants_file, 'w') as reduced_variants, open(filtered_variants_file, 'w') as filtered_variants:
        for line in f:
            variant = line.split()
            genotypes = [v.split(':')[0] for v in variant[3:]]
            depths = [v.split(':')[2] for v in variant[3:]]
            reference = variant[0][:10] + '8' + variant[0][11:]
            pos = variant[1]
            qual = variant[2]
            reduced_variants.write(reference + '\t' + pos + '\t' + qual + '\t' + '\t'.join(genotypes) + '\t' + '\t'.join(depths) + '\n')
            if float(qual) > 300 and any([test_depth(d) for d in depths]):
                filtered_variants.write(reference + '\t' + pos + '\t' + qual + '\t' + '\t'.join(genotypes) + '\t' + '\t'.join(depths) + '\n')