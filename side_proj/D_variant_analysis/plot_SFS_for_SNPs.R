library(MASS)
library(RColorBrewer)
library(AsexStats)

files <- paste0('/Volumes/dump/projects/timema/Timema_SNP_calling/data/', timemas$code, '_trinalge_SNP_filter_passed.tsv')
trinagles <- lapply(files, read.table, col.names = c('het', 'alleles', 'freq'))

get_decomp_tri <- function(tri){
    decomposed_trinagle <- matrix(0, nrow = 6, ncol = 10)
    for (i in 1:nrow(tri)){
        decomposed_trinagle[tri$het[i] + 1, tri$alleles[i]] = tri$freq[i]
    }
    decomposed_trinagle
}

decomposed_trinagles <- lapply(trinagles, get_decomp_tri)

plot_SFS <- function(decomposed_trinagle, main = 'Site freq spectra', col = "Grey"){
    site_freq <- colSums(decomposed_trinagle)

    barplot(site_freq, main = main, col = col)
}

for (i in 1:10) {
    pdf(paste0('figures/', timemas$code[i], '_SNP_SFS.pdf'))
    plot_SFS(decomposed_trinagles[[i]], timemas$code[i], ifelse(i %% 2 == 0, sex_red, asex_blue))
    dev.off()
}

#
# rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
#
# plot_SFS_vs_heterozygosity <- function(vcf_lines, main = 'Site freq spectra vs heterozygosity'){
#     site_freq <- apply(vcf_lines, 1, alt_alleles)
#     pop_heterozygosity <-  apply(vcf_lines, 1, gen_heterozygots)
#     cells <- table(paste(pop_heterozygosity, site_freq))
#     cell_coordinates <- lapply(strsplit(names(cells), ' '), as.numeric)
#
#     k <- kde2d(pop_heterozygosity, site_freq, n=40, h = 1)
#     image(k, col=rf(20)[1:15],
#           main = main,
#           xlab = 'number of heterozygous individuals',
#           ylab = 'allele frequency',
#           xlim = c(-0.1, 5.1), ylim = c(0.8, 10.2))
#     text(sapply(cell_coordinates, function(x){ if(x[1] == 0){ 0.1 } else { x[1]}} ), sapply(cell_coordinates, function(x){ x[2] } ), cells)
# }

