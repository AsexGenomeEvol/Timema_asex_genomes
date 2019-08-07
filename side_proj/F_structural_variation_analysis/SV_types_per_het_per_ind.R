## here the heterozygosity is population heterozygosity

library(AsexStats)
source('F_structural_variation_analysis/load_SV_calls.R')
source('F_structural_variation_analysis/filter_SV_calls.R')

# load SV calls
SV_call_manta_files <- paste0("data/", timemas$codes, "/variant_calls/", timemas$codes, "_survivor_manta_calls_union.vcf")
SV_calls <- load_SV_calls(SV_call_manta_files)

source('F_structural_variation_analysis/plot_SV_barplots.R')
by_type_tables <- get_by_type_SV_tables(SV_calls, filter_rare = T)
sv_genotype_tables <- lapply(SV_calls, get_sv_table, genotypes = T, T)
sv_genotype_tables <- lapply(sv_genotype_tables, function(SV_table) {     homoz <- rowSums(SV_table[,c('00','01','02','03','04','05')] == "1/1")
    het <- rowSums(SV_table[,c('00','01','02','03','04','05')] == "0/1")
    SV_table$homo <- homoz > 0 & het == 0
    SV_table$het <- het > 0 & homoz == 0
    SV_table$rare =  homoz + het == 1
    SV_table } )

# population_heteroz <- data.frame(sp = timemas$codes,
#            homo = round(sapply(sv_genotype_tables, function(x){ mean(x$homo) }), 4),
#            het = round(sapply(sv_genotype_tables, function(x){ mean(x$het) }), 4),
#            both = round(sapply(sv_genotype_tables, function(x){ 1 - mean(x$homo + x$het) }), 4))
# population_heteroz$rare_het <- sapply(sv_genotype_tables, function(SV_table){ sum(SV_table$rare & SV_table$het) / sum(SV_table$het) } )
# range(population_heteroz$rare_het[seq(2,10, by = 2)])

pdf('figures/SV_overview_zoomed_manta.pdf', width = 10, height = 8)
par(mfrow = c(4,1))
for ( type in c('DEL', 'DUP', 'INS', 'INV') ){
    ymax <- c(5000, 2000, 800, 800)[type == c('DEL', 'DUP', 'INS', 'INV')]
    main <- c('Deletions', 'Duplications', 'Insertions', 'Inversions')[type == c('DEL', 'DUP', 'INS', 'INV')]
    legend <- ifelse(type == 'DEL', T, F)
    plot_barplots(by_type_tables, type, main, ymax, legend)
}
dev.off()

pdf('figures/SV_overview_manta.pdf', width = 20, height = 16)
par(mfrow = c(4,1))
for ( type in c('DEL', 'DUP', 'INS', 'INV') ){
    main <- c('Deletions', 'Duplications', 'Insertions', 'Inversions')[type == c('DEL', 'DUP', 'INS', 'INV')]
    legend <- ifelse(type == 'DEL', T, F)
    plot_barplots(by_type_tables, type, "", NA, legend)
}
dev.off()

SV_call_lumpy_files <- paste0("data/", timemas$codes, "/variant_calls/", timemas$codes, "_survivor_lumpy_calls_union.vcf")
lumpy_SV_calls <- load_SV_calls(SV_call_lumpy_files)
lumpy_by_type_tables <- get_by_type_SV_tables(lumpy_SV_calls)

pdf('figures/SV_overview_lmupy_zoomed.pdf', width = 10, height = 6)
par(mfrow = c(3,1))
for ( type in c('DEL', 'DUP', 'INV') ){
    ymax <- c(5000, 2000, 800)[type == c('DEL', 'DUP', 'INV')]
    main <- c('Deletions', 'Duplications', 'Inversions')[type == c('DEL', 'DUP', 'INV')]
    legend <- ifelse(type == 'DEL', T, F)
    plot_barplots(lumpy_by_type_tables, type, main, ymax, legend)
}
dev.off()

pdf('figures/SV_overview_lmupy.pdf', width = 10, height = 6)
par(mfrow = c(3,1))
for ( type in c('DEL', 'DUP', 'INV') ){
    main <- c('Deletions', 'Duplications', 'Inversions')[type == c('DEL', 'DUP', 'INV')]
    legend <- ifelse(type == 'DEL', T, F)
    plot_barplots(lumpy_by_type_tables, type, main, NA, legend)
}
dev.off()

# Bah, this require a bit of thinking
# How many SVs I expect to find homo/heteroz in all given allelic frequency of the SV?
# get_AA <- function(p){
#     ((1 - p)^2 + 2 * p * (1 - p)) * ((p)^2)^5 +
#     (p^2 + 2 * p * (1 - p)) * ((1 - p)^2)^5
# }
# get_AB <- function(p){
#     2 * p * (1 - p) * ((p)^2)^5 + ...
# }
#
# ps <- seq(0.001, 0.999, length = 1000)
# plot(ps, get_AA(ps), type = 'l', ylim = c(0, 1))
# lines(ps, get_AB(ps), lty = 2)


#### final plot
source('F_structural_variation_analysis/plot_SV_barplots.R')

pdf('figures/SV_overview_manta_non_rare.pdf', width = 12, height = 10)
par(mfrow = c(3,1))
# 'DUP',
for ( type in c('DEL', 'INS', 'INV') ){
    plot_barplots(by_type_tables, type, "", NA, F)
}
dev.off()

pdf('figures/SV_invertions_verview_manta_non_rare_zoomed.pdf', width = 18, height = 10)
par(mfrow = c(1,1))
    plot_barplots(by_type_tables, "INV", "", 600, F)
dev.off()


pdf('figures/SV_merged_overview_manta.pdf', width = 20, height = 4)

source('F_structural_variation_analysis/plot_SV_barplots.R')
merged_types <- list()
merged_types[["MERGED"]] <- by_type_tables[[1]] + by_type_tables[[2]] + by_type_tables[[3]]
plot_barplots(merged_types, "MERGED", "", NA, F)

dev.off()