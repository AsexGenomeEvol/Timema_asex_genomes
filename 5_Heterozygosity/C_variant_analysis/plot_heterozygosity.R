library('RColorBrewer')
library('AsexStats')
# library('lme4')

heterozygosity_table_filename <- 'tables/heterozygosity_table.tsv'

all_genotypes <- paste0("0", 1:5, "_SNPs_het")
heterozygosity_table <- read.table(heterozygosity_table_filename, header = T, check.name = F)
rownames(heterozygosity_table) <- heterozygosity_table$sp

called_snps <- t(heterozygosity_table[,all_genotypes])

heterozygosities <- 100 * t(t(called_snps) / heterozygosity_table$callable_sites)

mins <- apply(heterozygosities, 2, min, na.rm = T)
maxes <- apply(heterozygosities, 2, max, na.rm = T)
bar_sizes <- colMeans(heterozygosities, na.rm = T)

library(reshape2)
individuals_het_tab <- melt(heterozygosities)
colnames(individuals_het_tab) <- c('ind', 'sp', 'het')

individuals_het_tab <- individuals_het_tab[!is.na(individuals_het_tab$het), ]
individuals_het_tab$pair <- substr(individuals_het_tab$sp, 1, 1)
individuals_het_tab$repr <- 'asex'

individuals_het_tab[individuals_het_tab$sp %in% timemas$codes[seq(2,10, by = 2)], 'repr'] <- 'sex'
# repr_mode_test <- glm(het ~ pair + sp + repr, data = individuals_het_tab)
# repr_mode_test <- lmer(het ~ pair + repr + (1| sp), data = individuals_het_tab, family = 'Guassian')
# repr_mode_test <- lme(het ~ pair + repr, random = ~1|sp, data=individuals_het_tab)

# anova(repr_mode_test)


### GenomeScope data

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

pdf('figures/heterozygosity_genomescope_SNPs_SVs.pdf', width = 10, height = 6)

par(mfrow = c(1, 2))

#### SNPs
ymax <- max(genomescope_est, bar_sizes, na.rm = T) + 0.05

locations <- barplot(genomescope_est, col = c(asex_blue, sex_red),
                     ylim = c(0, ymax), xaxt = "n", cex.axis = 1.3) # ylab = 'Number of called variants'
text(locations,
     par("usr")[3] - (0.03 * ymax), pos = 1, srt = 25,
     xpd = TRUE, labels = timemas$labels, cex = 0.8)
mtext('Nucleotide heterozygosity [%]', 2, padj = -3.2, cex = 1.3)

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
legend('topleft', pch = c(15, 20), col = c(asex_blue, sex_red), c('parthenogenic', 'sexual'), bty = 'n')

#### SVs

all_genotypes <- paste0("0", 1:5, "_SVs_het")
heterozygosity_table <- read.table(heterozygosity_table_filename, header = T, check.name = F)
rownames(heterozygosity_table) <- heterozygosity_table$sp

called_svs <- t(heterozygosity_table[,all_genotypes])
sv_heterozygosities <- t(t(called_svs) / heterozygosity_table$callable_sites)

mins <- apply(sv_heterozygosities, 2, min, na.rm = T)
maxes <- apply(sv_heterozygosities, 2, max, na.rm = T)
bar_sizes <- colMeans(sv_heterozygosities, na.rm = T)

ymax <- max(maxes, na.rm = T)

locations <- barplot(bar_sizes, col = c(asex_blue, sex_red),
                     ylim = c(0, ymax), xaxt = "n", cex.axis = 1.3) # ylab = 'Number of called variants'
text(locations,
     par("usr")[3] - (0.03 * ymax), pos = 1, srt = 25,
     xpd = TRUE, labels = timemas$labels, cex = 0.8)

mtext('SV heterozygosity', 2, padj = -3.2, cex = 1.3)

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

dev.off()
