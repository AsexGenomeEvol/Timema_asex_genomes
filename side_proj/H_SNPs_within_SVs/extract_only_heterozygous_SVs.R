########
library(AsexStats)

# load SV tables
source('F_structural_variation_analysis/load_SV_calls.R')
source('F_structural_variation_analysis/filter_SV_calls.R')

SV_call_manta_files <- paste0("data/", timemas$codes, "/variant_calls/", timemas$codes, "_survivor_manta_calls_union.vcf")
SV_calls <- load_SV_calls(SV_call_manta_files)

sv_genotype_tables <- lapply(SV_calls, get_sv_table, genotypes = T, T)

sv_tab <- sv_genotype_tables[[1]]

getNonRareHeterozygous <- function(sv_tab){
    het_rows <- apply(sv_tab[,4:9] == './.' | sv_tab[,4:9] == '0/1', 1, all) & (rowSums(sv_tab[,4:9] == '0/1') >= 2)
    sv_tab <- sv_tab[het_rows, c('scf', 'type', 'pos', 'end', 'len')]
    # kicking out translocations
    sv_tab <- sv_tab[sv_tab$type != "TRA", ]
    new_pos <- apply(sv_tab[, c('pos', 'end')], 1, min)
    new_end <- apply(sv_tab[, c('pos', 'end')], 1, max)
    sv_tab$pos <- new_pos
    sv_tab$end <- new_end
    sv_tab
}

het_SVs <- lapply(sv_genotype_tables, getNonRareHeterozygous)

SV_call_manta_heterozygous_variants_files <- paste0("data/", timemas$codes, "/variant_calls/", timemas$codes, "_survivor_manta_calls_union_non_rare_heterozygous_only.tsv")

for ( i in 1:10 ){
        write.table(het_SVs[[i]], SV_call_manta_heterozygous_variants_files[i], quote = F, sep = '\t', row.names = F)
}
