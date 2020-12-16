
library('AsexStats')

source('F_structural_variation_analysis/plot_SFS.R')


get_SFS <- function(sp){
    tab_filename <- paste0('data/SNP_calls/', sp, '_reduced_filtered_variants.tsv')
    variant_tab <- read.table(tab_filename, stringsAsFactors = F)
    colnames(variant_tab) <- c('scf', 'pos', 'qual', paste0('g', 1:5), paste0('d', 1:5), 'ref_scf', 'ref_pos', 'lg', 'lg_pos')

    genotypes <- paste0('g', 1:5)

    table(rowSums((1 * (variant_tab[,genotypes] == '0/1')) + (2 * (variant_tab[,genotypes] == '1/1'))))
}

Tms_SFS <- get_SFS('3_Tms')
Tce_SFS <- get_SFS('3_Tce')

par(mfrow = c(1, 2))

barplot(Tce_SFS, col = sex_red, ylim = c(0, max(Tce_SFS, Tms_SFS)))
legend('topright', title = 'T. cristinae', bty = 'n', legend = ' ')
barplot(Tms_SFS, col = asex_blue, ylim = c(0, max(Tce_SFS, Tms_SFS)))
legend('topright', title = "T. monikensis", bty = 'n', legend = ' ')

dev.off()


#####
# alt_alleles <- function(x){
#     sum(grepl("0/1", substr(x[10:15], 1, 3)) * 1) +
#     sum(grepl("1/1", substr(x[10:15], 1, 3)) * 2)
# }
#
# gen_heterozygots <- function(x){
#     sum(grepl("0/1", substr(x[10:15], 1, 3)))
# }
#
# plot_SFS <- function(vcf_lines, main = 'Site freq spectra', col = "Grey"){
#     site_freq <- sapply(vcf_lines, alt_alleles)
#
#     barplot(table(site_freq), main = main, col = col)
# }