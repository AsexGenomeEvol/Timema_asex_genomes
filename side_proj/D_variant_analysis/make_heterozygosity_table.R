library('RColorBrewer')
library('AsexStats')

heterozygosity_table_filename <- 'tables/heterozygosity_table.tsv'

if ( !file.exists(heterozygosity_table_filename)){

    heterozygosity_table <- data.frame(sp = timemas$code, row.names = timemas$code)

    all_genotypes <- c('g1', 'g2', 'g3', 'g4', 'g5')
    all_genotype_variants <- paste0(all_genotypes, '_variants')
    heterozygosity_table[, c(all_genotypes, all_genotype_variants, 'total_variants')] <- NA

    for(sp in timemas$code){
        tab_filename <- paste0('data/SNP_calls/', sp, '_reduced_filtered_variants.tsv')
        variant_tab <- read.table(tab_filename, stringsAsFactors = F)
        colnames(variant_tab) <- c('scf', 'pos', 'qual', paste0('g', 1:5), paste0('d', 1:5))

        if ( sp == '1_Tps'){
            genotypes <- c('g1', 'g2', 'g3')
        } else if ( sp == '4_Tte'){
            genotypes <- c('g2', 'g3', 'g4', 'g5')
        } else if ( sp == '2_Tsi'){
            genotypes <- c('g4', 'g5')
        } else {
            genotypes <- all_genotypes
        }
        genotype_variants <- paste0(genotypes, '_variants')

        heterozygous_variants <- colSums(variant_tab[,genotypes] == '0/1')
        called_variants <- colSums(variant_tab[,genotypes] != './.')

        heterozygosity_table[sp ,'total_variants'] <- nrow(variant_tab)
        heterozygosity_table[sp ,genotypes] <- heterozygous_variants
        heterozygosity_table[sp ,genotype_variants] <- called_variants
        # I could add asm stats to be able to divide by asm span
    }

    all_genotypes <- paste0("0", 1:5, "_SNPs_het")
    colnames(heterozygosity_table) <- c("sp", all_genotypes, paste0("0", 1:5, "_SNPs"), "total_sp_SNPs")

    # data from https://github.com/AsexGenomeEvol/Timema_SNP_calling/tree/master/F_snp_calling_round1#number-of-snp-and-monomorphic-positions-
    species <- c('1_Tdi', '1_Tps', '2_Tcm', '2_Tsi', '3_Tce', '3_Tms', '4_Tbi', '4_Tte', '5_Tge', '5_Tpa')
    callable_sites <- c(966024897, 874807942, 899565414, 896299726, 902715794, 991389399, 920445028, 980178552, 960991903, 705548015)
    heterozygosity_table[species, 'callable_sites'] <- callable_sites

    source('D_variant_analysis/load_SV_calls.R')

    # SVs_filt_relaxed_union_files <- paste0('data/manta_SV_calls/data/', timemas$code, '/SVs_filt_relaxed_union.vcf')
    SVs_filt_stringent_union_files <- paste0('data/manta_SV_calls/data/', timemas$code, '/SVs_filt_stringent_union.vcf')
    # SVs_filt_very_stringent_union_files <- paste0('data/manta_SV_calls/data/', timemas$code, '/SVs_filt_very_stringent_union.vcf')

    ### no filt
    SV_list <- load_SV_calls(SVs_filt_stringent_union_files)
    sv_genotype_tables <- lapply(SV_list, get_sv_table, genotypes = T)

    all_SV_genotypes <- c('01', '02', '03', '04', '05')
    for(i in 1:10){
        sp <- timemas$code[i]

        if ( sp == '1_Tps'){
            SV_genotypes <- c('01', '02', '03')
        } else if ( sp == '4_Tte'){
            SV_genotypes <- c('02', '03', '04', '05')
        } else if ( sp == '2_Tsi'){
            SV_genotypes <- c('04', '05')
        } else {
            SV_genotypes <- all_SV_genotypes
        }
        per_ind_heterozygous_SVs <- paste0(SV_genotypes, '_SVs_het')
        per_ind_SVs <- paste0(SV_genotypes, '_SVs')

        heterozygous_variants <- colSums(sv_genotype_tables[[i]][,SV_genotypes] == '0/1')
        called_variants <- colSums(sv_genotype_tables[[i]][,SV_genotypes] != './.')

        heterozygosity_table[sp , 'total_SVs'] <- nrow(sv_genotype_tables[[i]])
        heterozygosity_table[sp , per_ind_heterozygous_SVs] <- heterozygous_variants
        heterozygosity_table[sp , per_ind_SVs] <- called_variants
        # I could add asm stats to be able to divide by asm span
    }

    write.table(heterozygosity_table, 'tables/heterozygosity_table.tsv',  quote = F, sep = '\t', row.names = F)
} else {
    cat('Table exists already\n')
}