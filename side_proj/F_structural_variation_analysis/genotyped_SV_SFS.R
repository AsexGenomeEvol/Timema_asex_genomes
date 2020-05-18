# load genotype vcf file (either using the written functions, or even de novo)
# use F_structural_variation_analysis/plot_SFS.R functions for plotting SFS

library('AsexStats')

source('F_structural_variation_analysis/plot_SFS.R')

vcf_file <- 'data/genotyping/3_Tms_merged_calls_naive.vcf'
SVs <- strsplit(readLines(vcf_file), '\t')

SVs <- SVs[sapply(SVs, alt_alleles) > 0]

#### SV SFS

pdf('figures/genotyping/3_Tms_SV_genotyped_SFS.pdf')
    plot_SFS(SVs, "T. monikensis", asex_blue)
dev.off()

##### Triangle
source('F_structural_variation_analysis/plot_triangle.R')

site_freq <- sapply(SVs, alt_alleles)
pop_heterozygosity <- sapply(SVs, gen_heterozygots)

pdf('figures/genotyping/3_Tms_SV_genotyped_triangle.pdf')
    plot_trinagle(site_freq, pop_heterozygosity, asex_blue, "T. monikensis triangle")
dev.off()

for (type in c('INV', 'DUP', 'INS', 'DEL')){
    subset <- sapply(SVs, function(x){ grepl(type, x[3]) })
    site_freq <- sapply(SVs[subset], alt_alleles)
    pop_heterozygosity <- sapply(SVs[subset], gen_heterozygots)
    pdf(paste0('figures/genotyping/3_Tms_SV_genotyped_triangle_', type, '.pdf'))
        plot_trinagle(site_freq, pop_heterozygosity, asex_blue, paste("T. monikensis", type))
    dev.off()
}

# table of types
table(sapply(SVs, function(x){ substr(x[3], 6, 8)} ))
#  DEL  DUP  INS  INV
# 1463   57  510   12
# Oh, my. There you go inversions!

inversions <- SVs[sapply(SVs, function(x){ substr(x[3], 6, 8) == 'INV'} )]
# look at the overlaps of those and SNPs, there are so few of them that it should be very easy
# perhaps look at BAM files at the location of the inversions in genome browser


text <- sapply(inversions, function(x) {paste(x, collapse = '\t')} )
writeLines(text, con = 'data/genotyping/3_Tms_merged_calls_naive_INV.vcf', sep = "\n", useBytes = FALSE)

# > SO F* up
# solid ones:
# 3_Tms_b3v08_scaf002312: 21725 - 23286; 00
# 3_Tms_b3v08_scaf000092: 80652 - 80811; 00 01 02 05

