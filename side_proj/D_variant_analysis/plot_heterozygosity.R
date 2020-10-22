library('RColorBrewer')
library('AsexStats')

heterozygosity_table_filename <- 'tables/heterozygosity_table.tsv'

all_genotypes <- paste0("0", 1:5, "_SNPs_het")
heterozygosity_table <- read.table(heterozygosity_table_filename, header = T, check.name = F)
rownames(heterozygosity_table) <- heterozygosity_table$sp

called_snps <- t(heterozygosity_table[,all_genotypes])

heterozygosities <- 100 * t(t(called_snps) / heterozygosity_table$callable_sites)

mins <- apply(heterozygosities, 2, min, na.rm = T)
maxes <- apply(heterozygosities, 2, max, na.rm = T)
bar_sizes <- colMeans(heterozygosities, na.rm = T)

ymax <- max(maxes)

### plain
# pdf('figures/heterozygosity_reseq_data.pdf', width = 12, height = 6)
#
# locations <- barplot(bar_sizes, col = c(asex_blue, sex_red),
#                      ylim = c(0, ymax), xaxt = "n", cex.axis = 1.3) # ylab = 'Number of called variants'
# text(locations,
#      par("usr")[3] - 50, pos = 1,
#      xpd = TRUE, labels = timemas$labels, cex = 1)
# mtext('percentage of heterozygous callable sites', 2, padj = -3.2, cex = 1.3)
#
# w <- 0.1
# for ( i in 1:length(mins) ){
#     x <- locations[i]
#     y_min <- mins[i]
#     y_max <- maxes[i]
#     lines(c(x, x), c(y_min, y_max), xpd=T, lwd = 1.5)
#     lines(c(x - w, x + w), c(y_min, y_min), xpd=T, lwd = 1.5)
#     lines(c(x - w, x + w), c(y_max, y_max), xpd=T, lwd = 1.5)
# }
#
# dev.off()

### with GenomeScope
pdf('figures/heterozygosity_reseq_data_genomescope.pdf', width = 12, height = 6)

source('../timema_assembly/B2_genome_profiling/extract_heterozygosity_genomescope.R')
genomescope_files <- sapply(timemas$codes, function(x){(dir(paste0('../timema_assembly/data/',x,'/genomescope'), pattern = 'summary', full.names = T))})
genomescope_est_min <- sapply(genomescope_files, extract_heterozygosity, min = T, rounding = 6)
genomescope_est_max <- sapply(genomescope_files, extract_heterozygosity, min = F, rounding = 6)
genomescope_est <- (genomescope_est_min + genomescope_est_max) / 2
# > genomescope_est
#     4_Tte     4_Tbi     2_Tsi     2_Tcm     1_Tdi     1_Tps     3_Tms     3_Tce
# 0.0877690 0.3624495 0.0823455 0.6148230 0.1074665 0.4232380 0.0859420 0.9479600
#     5_Tge     5_Tpa
# 0.1029225 2.1619950
genomescope_est[seq(1,10, by = 2)] <- NA

ymax <- max(genomescope_est, bar_sizes, na.rm = T) + 0.05

locations <- barplot(genomescope_est, col = c(asex_blue, sex_red),
                     ylim = c(0, ymax), xaxt = "n", cex.axis = 1.3) # ylab = 'Number of called variants'
text(locations,
     par("usr")[3], pos = 1,
     xpd = TRUE, labels = timemas$labels, cex = 0.8)
mtext('Heterozygosity [%]', 2, padj = -3.2, cex = 1.3)

w <- 0.1
for ( i in 1:length(mins) ){
    x <- locations[i]
    y_min <- mins[i]
    y_max <- maxes[i]
    lines(c(x, x), c(y_min, y_max), xpd=T, lwd = 1.5)
    lines(c(x - w, x + w), c(y_min, y_min), xpd=T, lwd = 1.5)
    lines(c(x - w, x + w), c(y_max, y_max), xpd=T, lwd = 1.5)
}
points(bar_sizes ~ locations, pch = 20, cex = 0.5)
legend('topleft', pch = c(15, 20), col = c(sex_red, 1), c('GenomeScope', 'SNPs'), bty = 'n')

dev.off()

### SV heterozygosity
# source('D_variant_analysis/load_SV_calls.R')
#
# SVs_no_filt_union_files <- paste0('data/manta_SV_calls/data/', timemas$code, '/SVs_no_filt_union.vcf')
# SVs_filt_relaxed_union_files <- paste0('data/manta_SV_calls/data/', timemas$code, '/SVs_filt_relaxed_union.vcf')
# SVs_filt_stringent_union_files <- paste0('data/manta_SV_calls/data/', timemas$code, '/SVs_filt_stringent_union.vcf')
# SVs_filt_very_stringent_union_files <- paste0('data/manta_SV_calls/data/', timemas$code, '/SVs_filt_very_stringent_union.vcf')
#
# ### no filt
# SV_list <- load_SV_calls(SVs_filt_relaxed_union_files)
# sv_genotype_tables <- lapply(SV_list, get_sv_table, genotypes = T)
#
# all_SV_genotypes <- c('01', '02', '03', '04', '05')
# for(i in 1:10){
#     sp <- timemas$code[i]
#
#     if ( sp == '1_Tps'){
#         SV_genotypes <- c('01', '02', '03')
#     } else if ( sp == '4_Tte'){
#         SV_genotypes <- c('02', '03', '04', '05')
#     } else if ( sp == '2_Tsi'){
#         SV_genotypes <- c('04', '05')
#     } else {
#         SV_genotypes <- all_SV_genotypes
#     }
#     per_ind_heterozygous_SVs <- paste0(SV_genotypes, '_SVs_het')
#     per_ind_SVs <- paste0(SV_genotypes, '_SVs')
#
#     heterozygous_variants <- colSums(sv_genotype_tables[[i]][,SV_genotypes] == '0/1')
#     called_variants <- colSums(sv_genotype_tables[[i]][,SV_genotypes] != './.')
#
#     heterozygosity_table[sp , 'total_SVs'] <- nrow(sv_genotype_tables[[i]])
#     heterozygosity_table[sp , per_ind_heterozygous_SVs] <- heterozygous_variants
#     heterozygosity_table[sp , per_ind_SVs] <- called_variants
#     # I could add asm stats to be able to divide by asm span
# }
#
# per_ind_heterozygous_SVs <- paste0(all_SV_genotypes, '_SVs_het')
# per_ind_SVs <- paste0(all_SV_genotypes, '_SVs')
#
# called_snps <- t(heterozygosity_table[,per_ind_heterozygous_SVs])
#
# heterozygosities <- 100 * t(t(called_snps) / heterozygosity_table$callable_sites)
#
# mins <- apply(heterozygosities, 2, min, na.rm = T)
# maxes <- apply(heterozygosities, 2, max, na.rm = T)
# bar_sizes <- colMeans(heterozygosities, na.rm = T)
#
# ymax <- max(maxes)
#
# ### plain
# pdf('figures/heterozygosity_SV_reseq_data_FILT_RELAXED.pdf', width = 12, height = 6)
#
# locations <- barplot(bar_sizes, col = c(asex_blue, sex_red),
#                      ylim = c(0, ymax), xaxt = "n", cex.axis = 1.3, main = 'SVs heterozygosity - relaxed') # ylab = 'Number of called variants'
# text(locations,
#      par("usr")[3] - 50, pos = 1,
#      xpd = TRUE, labels = timemas$labels, cex = 1)
# mtext('percentage of heterozygous callable sites', 2, padj = -3.2, cex = 1.3)
#
# w <- 0.1
# for ( i in 1:length(mins) ){
#     x <- locations[i]
#     y_min <- mins[i]
#     y_max <- maxes[i]
#     lines(c(x, x), c(y_min, y_max), xpd=T, lwd = 1.5)
#     lines(c(x - w, x + w), c(y_min, y_min), xpd=T, lwd = 1.5)
#     lines(c(x - w, x + w), c(y_max, y_max), xpd=T, lwd = 1.5)
# }
#
# dev.off()