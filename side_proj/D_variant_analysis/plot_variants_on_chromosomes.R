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
        pal <- brewer.pal(5, "YlGnBu")[c(3,5)]
    } else {
        pal <- brewer.pal(5, "YlOrRd")[c(3,5)]
    }

    plot(NULL, xlim = c(1, nrow(variant_density_table)), ylim = c(0, 1), pch = 20,
         main = sp, xaxt = "n", yaxt = "n", bty = 'n', xlab = '', ylab = '', cex.axis = 1.4, cex.main = 1.6)
         # xlab = 'linage group [ Mbp ]', ylab = '# found SVs'
    xtick <- chromosomes$adjustments / window
    axis(side = 1, at = xtick, labels = FALSE)
    text(x = (chromosomes$adjustments[1:12] + chromosomes$adjustments[2:13]) / (2 * window),  par("usr")[3],
         labels = chromosomes$chr[1:12], pos = 1, xpd = TRUE, cex = 1.3)
    for (i in seq(1, 12, by = 2)) {
        rect(chromosomes$adjustments[i] / window, -2, chromosomes$adjustments[i + 1] / window, 1, col = 'lightgrey', border = F)
    }
    if ( sex == 'asex'){
        axis(side = 2, at = seq(0,1, by = 0.2), labels = seq(0, SNP_ylim, length = 6))
        mtext('SNPs', 2, line = 2.5)
    } else {
        axis(side = 4, at = seq(0,1, by = 0.2), labels = seq(0, SV_ylim, length = 6))
        mtext('SVs', 4, line = 2.5)
    }

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
    legend('topright', bty = 'n', c('SNPs', 'SVs'), pch = 20, col = pal[c(1,2)])
}

##########

source('D_variant_analysis/load_chromosomes.R')

pdf('figures/anchored_SNPs_and_SVs_overlayed.pdf', width = 10, height = 6)
par(mfrow = c(5, 2))
for(sp in timemas$codes){

    variant_density_table_file <- paste0("tables/", sp, "_variants_on_chromosomes_w", window, ".tsv")

    if ( file.exists(variant_density_table_file) ){
        variant_density_table <- read.table(variant_density_table_file, header = T, sep = '\t', stringsAsFactors = F)
    } else {
        print('please, run D_variant_analysis/make_variant_density_table.R first')
    }

    if ( sp %in% timemas$codes[seq(1, 10, by=2)]){
        par(mar=c(2,4,2,0))
        # 'c(bottom, left, top, right)'
        plot_SVs_on_LGs(F, T, 'asex', 28000, 100)
    } else {
        par(mar=c(2,0,2,4))
        plot_SVs_on_LGs(F, T, 'sex', 28000, 100)
    }

}

dev.off()