library(MASS)
library(RColorBrewer)

alt_alleles <- function(x){
    sum(grepl("0/1", substr(x[10:15], 1, 3)) * 1) +
    sum(grepl("1/1", substr(x[10:15], 1, 3)) * 2)
}

gen_heterozygots <- function(x){
    sum(grepl("0/1", substr(x[10:15], 1, 3)))
}

plot_SFS <- function(vcf_lines, main = 'Site freq spectra', col = "Grey"){
    site_freq <- sapply(vcf_lines, alt_alleles)

    barplot(table(site_freq), main = main, col = col)
}
