library(AsexStats)
library(RColorBrewer)

# asexuals <- seq(1, 10, by=2)
# [asexuals]
to_analyse <- timemas$codes


### LOADING SNPs
hetero_SNP_files <- paste0('Timema_SNP_calling/data/', to_analyse, '_heterozygous_SNP_filter_passed.tsv')
homo_SNP_files <- paste0('Timema_SNP_calling/data/', to_analyse, '_homozygous_SNP_filter_passed.tsv')

heteroz_variants <- lapply(hetero_SNP_files, read.table, col.names=c('scf', 'pos'), stringsAsFactors = F)
homo_variants <- lapply(homo_SNP_files, read.table, col.names=c('scf', 'pos'), stringsAsFactors = F)

# loading SVs
SV_call_manta_heterozygous_variants_files <- paste0("data/", to_analyse, "/variant_calls/", to_analyse, "_survivor_manta_calls_union_non_rare_heterozygous_only.tsv")
hetero_SVs <- lapply(SV_call_manta_heterozygous_variants_files, read.table, header = T, stringsAsFactors = F)

# one_snp <- heteroz_variants[[1]]
# one_sv <- hetero_SVs[[1]]

### get simple number of scaffolds that has both heterozygous SVs and SNPs
get_total_scfs_and_shared <- function(one_snp, one_sv){
    tab <- table(table(c(unique(one_snp$scf), unique(one_sv$scf))))
    c(sum(tab), tab[2])
}

sapply(1:5, function(i){ get_total_scfs_and_shared(heteroz_variants[[i]], hetero_SVs[[i]]) } )
# heterozygous SVs or SNPs are found on 3746 scaffolds, of which both SVs and SNPs contain 556
# similar ratios for others too
#   [,1] [,2] [,3] [,4] [,5]
#   3746 4355 4591 2981 2876 (total)
# 2  556  513  613  575  370 (both SNPs and SVs)

### extract those that contain both

get_list_of_SNPs_on_SVs <- function(one_snp, one_sv, margin = 500, INV_only = F){
    tab_of_scfs <- table(c(unique(one_snp$scf), unique(one_sv$scf)))
    both_SNPs_and_SVs <- names(tab_of_scfs)[tab_of_scfs == 2]

    ### let's go through these and find if there are any overlaps
    SNPs_on_SVs <- data.frame()
    SNPs_on_INV <- data.frame()

    for (scf in both_SNPs_and_SVs){
        scf_snps <- one_snp[one_snp$scf == scf,]
        scf_SVs <- one_sv[one_sv$scf == scf,]
        scf_INV <- one_sv[one_sv$scf == scf & one_sv$type == 'INV',]
        for (snp_i in 1:nrow(scf_snps)) {
            snp <- scf_snps[snp_i, ]
            if ( any(snp$pos > (scf_SVs$pos - margin) & snp$pos < (scf_SVs$end + margin)) ){
                SNPs_on_SVs <- rbind(SNPs_on_SVs, snp)
            }
            if ( any(snp$pos > (scf_INV$pos - margin) & snp$pos < (scf_INV$end + margin)) ){
                SNPs_on_INV <- rbind(SNPs_on_INV, snp)
            }
        }
    }
    if ( INV_only ){
        SNPs_on_INV
    } else {
        SNPs_on_SVs
    }
}

SNPs_directly_on_INVs <- lapply(1:5, function(i){ get_list_of_SNPs_on_SVs(heteroz_variants[[i]], hetero_SVs[[i]], 0, T) } )
SNPs_around_INVs <- lapply(1:5, function(i){ get_list_of_SNPs_on_SVs(heteroz_variants[[i]], hetero_SVs[[i]], 10000, T) } )
SNPs_around_SVs <- lapply(1:5, function(i){ get_list_of_SNPs_on_SVs(heteroz_variants[[i]], hetero_SVs[[i]], 10000, F) } )

summary_tab <- data.frame(
    sp = to_analyse,
    total_heterozygous_SNPs = sapply(heteroz_variants, nrow),
    total_heterozygous_SVs = sapply(hetero_SVs, nrow),
    SNPs_directly_on_INV = sapply(SNPs_directly_on_INVs, nrow),
    SNPs_around_INV = sapply(SNPs_around_INVs, nrow),
    SNPs_around_SVs = sapply(SNPs_around_SVs, nrow))

write.table(summary_tab, 'tables/SNP_SV_overlap_asexuals.tsv')
