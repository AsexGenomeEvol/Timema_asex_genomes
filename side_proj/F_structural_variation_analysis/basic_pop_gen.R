library(AsexStats)
source('F_structural_variation_analysis/plot_SFS.R')
source('F_structural_variation_analysis/load_SV_calls.R')

###############
# MANTA CALLS #
###############

SV_call_manta_files <- paste0("data/", timemas$codes, "/variant_calls/", timemas$codes, "_survivor_manta_calls_union.vcf")
SV_calls <- load_SV_calls(SV_call_manta_files)

# SFSs
pdf('figures/SFS_manta.pdf', width = 12, height = 8)
layout(mat = matrix(1:10, nrow = 2, ncol = 5))
for( i in 1:10){
    plot_SFS(SV_calls[[i]], timemas$labels[i], timemas$cols[i])
}
dev.off()

pdf('figures/SFS_vs_heterozygosity_manta.pdf', width = 30, height = 12)
layout(mat = matrix(1:10, nrow = 2, ncol = 5))
for( i in 1:10){
    plot_SFS_vs_heterozygosity(SV_calls[[i]], timemas$labels[i])
}
dev.off()

### Tms individual plots

i <- 8
Tce_site_freq <- sapply(SV_calls[[i]], alt_alleles)
Tce_pop_heterozygosity <- sapply(SV_calls[[i]], gen_heterozygots)

pdf('figures/SFS_Tce_manta.pdf', width = 9, height = 4)
    barplot(table(Tce_site_freq), col = sex_red, cex.axis = 2, cex.names = 1.75)
dev.off()

cells <- table(paste(Tce_pop_heterozygosity, Tce_site_freq))
cell_coordinates <- lapply(strsplit(names(cells), ' '), as.numeric)
cr = colorRampPalette(c(sex_red, sex_red, sex_red, sex_red, 'white'))
k <- kde2d(Tce_site_freq, Tce_pop_heterozygosity, n = 100)

pdf('figures/SFS_vs_heterozygosity_Tce_manta.pdf', width = 6, height = 6)
    image(k, col = rev(cr(1000)),
          # ylab = 'number of heterozygous individuals',
          # xlab = 'allele frequency',
          ylim = c(-0.2, 6.2), xlim = c(0.6, 12.4),
          cex.axis = 1.4, bty = 'n')
    text(sapply(cell_coordinates, function(x){ if( x[2] == 1) { 1.2 } else { x[2] } } ),
         sapply(cell_coordinates, function(x){ x[1] } ),
         cells)
dev.off()

### Tce individual plots

i <- 7
Tms_site_freq <- sapply(SV_calls[[i]], alt_alleles)
Tms_pop_heterozygosity <- sapply(SV_calls[[i]], gen_heterozygots)

pdf('figures/SFS_Tms_manta.pdf', width = 9, height = 4)
    barplot(table(Tms_site_freq), col = asex_blue, cex.axis = 2, cex.names = 1.75)
dev.off()

cells <- table(paste(Tms_pop_heterozygosity, Tms_site_freq))
cell_coordinates <- lapply(strsplit(names(cells), ' '), as.numeric)
cr = colorRampPalette(c(asex_blue, asex_blue, asex_blue, asex_blue, 'white'))
k <- kde2d(Tms_site_freq, Tms_pop_heterozygosity, n = 100)

pdf('figures/SFS_vs_heterozygosity_Tms_manta.pdf', width = 6, height = 6)
    image(k, col = rev(cr(1000)),
          ylim = c(-0.2, 6.2), xlim = c(0.6, 12.4),
          cex.axis = 1.4, bty = 'n')
    text(sapply(cell_coordinates, function(x){ if( x[2] == 1) { 1.2 } else { x[2] } } ),
         sapply(cell_coordinates, function(x){ x[1] } ),
         cells)
dev.off()

####
# some more exploration

sv_genotype_tables <- lapply(SV_calls, get_sv_table, genotypes = T)

SVs <- sv_genotype_tables[[1]]

sv_genotype_tables <- lapply(sv_genotype_tables, function(SVs){
homoz <- rowSums(SVs[,c('00','01','02','03','04','05')] == "1/1")
het <- rowSums(SVs[,c('00','01','02','03','04','05')] == "0/1")
SVs$num_homoz <- homoz
SVs$num_het <- het
SVs$homo <- homoz > 0 & het == 0
SVs$het <- het > 0 & homoz == 0
SVs$both <- het > 0 & homoz > 0
SVs })

homo_non_rare <- sapply(sv_genotype_tables, function(x) { sum(x$homo & x$num_homoz > 1) / sum(x$homo) } )[seq(1,10,by=2)]
het_non_rare  <- sapply(sv_genotype_tables, function(x) { sum(x$het & x$num_het > 1) / sum(x$het) } )[seq(1,10,by=2)]
plot(homo_non_rare ~ het_non_rare, pch = 20, col = 1, bty = 'n', cex = 2, xlim = c(0.2, 0.7), ylim = c(0.2, 0.7),
xlab = 'Heterozygous in more than one ind [%]', ylab = 'Homozygous in more than one ind [%]')
lines(c(0, 1), c(0, 1))
points(homo_non_rare ~ het_non_rare, pch = 20, col = asex_blue, cex = 1.7)
text(het_non_rare, homo_non_rare, timemas$names[seq(1,10,by=2)], pos = 2)

# among asexuals
range(sapply(sv_genotype_tables, nrow)[seq(1,10, by = 2)])
# among sexuals
range(sapply(sv_genotype_tables, nrow)[seq(2,10, by = 2)])

sapply(sv_genotype_tables, function(SVs) { sum(SVs$homo) })
sapply(sv_genotype_tables, function(SVs) { sum(SVs$het) })
sapply(sv_genotype_tables, function(SVs) { sum(SVs$both) })

################
# SMOOVE CALLS #
################

SV_call_lumpy_files <- paste0("data/", timemas$codes, "/variant_calls/", timemas$codes, "_survivor_lumpy_calls_union.vcf")
SV_calls <- load_SV_calls(SV_call_lumpy_files)
# species list of SVs lists

# SFSs
pdf('figures/SFS_smove.pdf', width = 12, height = 8)
layout(mat = matrix(1:10, nrow = 2, ncol = 5))
for( i in 1:10){
    plot_SFS(SV_calls[[i]], timemas$labels[i], timemas$cols[i])
}
dev.off()

pdf('figures/SFS_vs_heterozygosity_smove.pdf', width = 30, height = 12)
layout(mat = matrix(1:10, nrow = 2, ncol = 5))
for( i in 1:10){
    plot_SFS_vs_heterozygosity(SV_calls[[i]], timemas$labels[i])
}
dev.off()
