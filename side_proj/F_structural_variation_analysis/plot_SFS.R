library(MASS)
library(RColorBrewer)

alt_alleles <- function(x){
    sum(grepl("0/1", x[10:15]) * 1) +
    sum(grepl("1/1", x[10:15]) * 2)
}

gen_heterozygots <- function(x){
    sum(grepl("0/1", x[10:15]))
}

plot_SFS <- function(vcf_lines, main = 'Site freq spectra', col = "Grey"){
    site_freq <- sapply(vcf_lines, alt_alleles)

    barplot(table(site_freq), main = main, col = col)
}

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))

plot_SFS_vs_heterozygosity <- function(vcf_lines, main = 'Site freq spectra vs heterozygosity'){
    site_freq <- sapply(vcf_lines, alt_alleles)
    pop_heterozygosity <- sapply(vcf_lines, gen_heterozygots)
    cells <- table(paste(pop_heterozygosity, site_freq))
    cell_coordinates <- lapply(strsplit(names(cells), ' '), as.numeric)

    k <- kde2d(pop_heterozygosity, site_freq, n=40)
    image(k, col=rf(20)[1:15],
          main = main,
          xlab = 'number of heterozygous individuals',
          ylab = 'allele frequency',
          xlim = c(-0.1, 6.1), ylim = c(0.8, 12.2))
    text(sapply(cell_coordinates, function(x){ if(x[1] == 0){ 0.1 } else { x[1]}} ), sapply(cell_coordinates, function(x){ x[2] } ), cells)
}

