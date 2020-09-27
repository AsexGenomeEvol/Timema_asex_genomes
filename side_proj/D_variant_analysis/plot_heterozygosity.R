library('RColorBrewer')
library('AsexStats')


heterozygosity_table <- data.frame(sp = timemas$code, row.names = timemas$code)

all_genotypes <- c('g1', 'g2', 'g3', 'g4', 'g5')
all_genotype_variants <- paste0(all_genotypes, '_variants')
heterozygosity_table[, c(all_genotypes, all_genotype_variants, 'total_variants')] <- NA

for(sp in timemas$code){
    tab_filename <- paste0('data/SNP_calls/', sp, '_reduced_filtered_variants.tsv')
    variant_tab <- read.table(tab_filename, stringsAsFactors = F)
    colnames(variant_tab) <- c('scf', 'pos', 'qual', paste0('g', 1:5), paste0('d', 1:5))

    if ( sp == '1_Tps'){
        genotypes <- c('g1', 'g2', 'g3')
    } else if ( sp == '4_Tte'){
        genotypes <- c('g2', 'g3', 'g4', 'g5')
    } else if ( sp == '2_Tsi'){
        genotypes <- c('g4', 'g5')
    } else {
        genotypes <- all_genotypes
    }
    genotype_variants <- paste0(genotypes, '_variants')

    heterozygous_variants <- colSums(variant_tab[,genotypes] == '0/1')
    called_variants <- colSums(variant_tab[,genotypes] != './.')

    heterozygosity_table[sp ,'total_variants'] <- nrow(variant_tab)
    heterozygosity_table[sp ,genotypes] <- heterozygous_variants
    heterozygosity_table[sp ,genotype_variants] <- called_variants
    # I could add asm stats to be able to divide by asm span
}

heterozygosities <- t(heterozygosity_table[,all_genotypes])

mins <- apply(heterozygosities, 2, min, na.rm = T)
maxes <- apply(heterozygosities, 2, max, na.rm = T)
bar_sizes <- colMeans(heterozygosities, na.rm = T)

ymax <- max(maxes)


pdf('figures/heterozygosity_reseq_data.pdf', width = 12, height = 6)

locations <- barplot(bar_sizes, col = c(asex_blue, sex_red),
                     ylim = c(0, ymax), xaxt = "n", cex.axis = 1.3) # ylab = 'Number of called variants'
text(locations,
     par("usr")[3] - 50, pos = 1,
     xpd = TRUE, labels = timemas$labels, cex = 1,
     ylab = 'number of heterozygous loci')
mtext("number of called heterozygous alleles", 2, padj = -3.2, cex = 1.3)

w <- 0.1
for ( i in 1:length(mins) ){
    x <- locations[i]
    y_min <- mins[i]
    y_max <- maxes[i]
    lines(c(x, x), c(y_min, y_max), xpd=T, lwd = 1.5)
    lines(c(x - w, x + w), c(y_min, y_min), xpd=T, lwd = 1.5)
    lines(c(x - w, x + w), c(y_max, y_max), xpd=T, lwd = 1.5)
}

dev.off()