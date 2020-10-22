library(AsexStats)
library(RColorBrewer)

###############################
##### plotting parameters #####
###############################

window = 1e6
gap_beween_chromosomes = 3

###############################

plot_SVs_on_LGs <- function(exclude_zero = T, lines = F, sex = 'asex', SNP_ylim, SV_ylim){

    if ( sex == 'asex' ){
        pair_index <- (which(sp == timemas$codes) + 1) / 2
        pal <- brewer.pal(5, "YlGnBu")[c(3,5)]
        plot(NULL, xlim = c(1, nrow(variant_density_table)), ylim = c(0, 1), pch = 20,
             main = timema_pairs$labels[pair_index],
             xaxt = "n", yaxt = "n", bty = 'n', xlab = '', ylab = '', cex.axis = 1.4, cex.main = 1.6)
    } else {
        pal <- brewer.pal(5, "YlOrRd")[c(3,5)]
        plot(NULL, xlim = c(1, nrow(variant_density_table)), ylim = c(0, 1), pch = 20,
             xaxt = "n", yaxt = "n", bty = 'n', xlab = '', ylab = '', cex.axis = 1.4, cex.main = 1.6)
    }

     # main = timemas$labels[timemas$codes == sp],
     # xlab = 'linage group [ Mbp ]', ylab = '# found SVs'
    xtick <- chromosomes$adjustments / window
    # axis(side = 1, at = xtick, labels = FALSE)
    if ( sp == '5_Tpa'){
        text(x = (chromosomes$adjustments[1:12] + chromosomes$adjustments[2:13]) / (2 * window),  par("usr")[3],
             labels = 1:12, pos = 1, xpd = TRUE, cex = 1, xpd=NA)
        text(x = (chromosomes$adjustments[11] + chromosomes$adjustments[12]) / (2 * window) + 50,  par("usr")[3],
             labels = 12, pos = 1, xpd = TRUE, cex = 1, xpd=NA)
    }

    for (i in seq(1, 12, by = 2)) {
        rect(chromosomes$adjustments[i] / window, -2, chromosomes$adjustments[i + 1] / window, 1, col = 'lightgrey', border = F)
    }
    # if ( sex == 'asex'){
        axis(side = 2, at = seq(0,1, by = 0.2), labels = seq(0, SNP_ylim, length = 6))
        mtext('SNPs', 2, line = 2.5)
    # } else {
        axis(side = 4, at = seq(0,1, by = 0.2), labels = seq(0, SV_ylim, length = 6))
        mtext('SVs', 4, line = 2.5)
    # }

    if ( exclude_zero ){
        variant_density_table[variant_density_table$SNPs == 0, 'SNPs'] <- NA
        variant_density_table[variant_density_table$SVs == 0, 'SVs'] <- NA
    }

    if ( lines ){
        lines(variant_density_table$SVs / SV_ylim, col = pal[2], lwd = 1.6)
        lines(variant_density_table$SNPs / SNP_ylim, col = pal[1], lwd = 1.6)
    } else {
        points(variant_density_table$SVs / SV_ylim, pch = 20, col = pal[2])
        points(variant_density_table$SNPs / SNP_ylim, pch = 20, col = pal[1])
    }
    # legend('topright', bty = 'n', c('SNPs', 'SVs'), pch = 20, col = pal[c(1,2)])
}

##########

source('D_variant_analysis/load_chromosomes.R')

variant_density_table_files <- paste0("tables/", timemas$codes, "_variants_on_chromosomes_w", window, ".tsv")
variant_density_tables <- lapply(variant_density_table_files, read.table, header = T, sep = '\t', stringsAsFactors = F)

SNPs_max <- max(sapply(variant_density_tables, function(x){ x[, 'SNPs'] } ))
SVs_max <- max(sapply(variant_density_tables, function(x){ x[, 'SVs'] } ))
print(SNPs_max)
print(SVs_max)

pdf('figures/anchored_SNPs_and_SVs_overlayed.pdf', width = 6, height = 10)
par(mfrow = c(10, 1), oma = c(4,0,2,0))
for(i in 1:10){
    sp = timemas$codes[i]
    variant_density_table <- variant_density_tables[[i]]

    if ( sp %in% timemas$codes[seq(1, 10, by=2)]){
        # 'c(bottom, left, top, right)'
        par(mar=c(0,4,1.25,4))
        plot_SVs_on_LGs(F, T, 'asex', 25500, 100)
    } else {
        par(mar=c(1.25,4,0,4))
        plot_SVs_on_LGs(F, T, 'sex', 25500, 100)
    }
    mtext("Linkage Group", side=1, line=2, cex=1.2, outer=TRUE)
}

dev.off()