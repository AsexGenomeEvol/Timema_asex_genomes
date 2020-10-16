library(AsexStats)
library(RColorBrewer)

###############################
##### plotting parameters #####
###############################

window = 1e6
gap_beween_chromosomes = 3

###############################

plot_SVs_on_LGs <- function(exclude_zero = T, lines = F, pal = 'asex', ylim = NA){

    # BrBG
    if ( pal == 'asex' ){
        pal <- brewer.pal(5, "YlGnBu")[c(3,5)]
    } else {
        pal <- brewer.pal(5, "YlOrRd")[c(3,5)]
    }
    # pal <- addalpha(pal, alpha)
    if ( any(is.na(ylim)) ){
        ylim = range(variant_density_table$variants, na.rm = T)
    }
    # print(ylim)

    plot(NULL, xlim = c(1, nrow(variant_density_table)), ylim = ylim, pch = 20,
         main = sp, xaxt = "n", bty = 'n', xlab = '', ylab = '', cex.axis = 1.4, cex.main = 1.6)
         # xlab = 'linage group [ Mbp ]', ylab = '# found SVs'
    xtick <- chromosomes$adjustments / window
    axis(side = 1, at = xtick, labels = FALSE)
    text(x = (chromosomes$adjustments[1:12] + chromosomes$adjustments[2:13]) / (2 * window),  par("usr")[3],
         labels = chromosomes$chr[1:12], pos = 1, xpd = TRUE, cex = 1.3)
    for (i in seq(1, 12, by = 2)) {
        rect(chromosomes$adjustments[i] / window, -2, chromosomes$adjustments[i + 1] / window, ylim[2], col = 'lightgrey', border = F)
    }

    if ( exclude_zero ){
        variant_density_table[variant_density_table$variants == 0, 'variants'] <- NA
        variant_density_table[variant_density_table$SVs == 0, 'common'] <- NA
    }

    if ( lines ){
        lines(variant_density_table$SVs / max(variant_density_table$SVs), col = pal[2], lwd = 1.6)
        lines(variant_density_table$variants / max(variant_density_table$variants), col = pal[1], lwd = 1.6)
    } else {
        points(variant_density_table$variants, pch = 20, col = pal[1])
        points(variant_density_table$SVs, pch = 20, col = pal[2])
    }
    legend('topright', bty = 'n', c('SNPs', 'SVs'), pch = 20, col = pal[c(1,2)])
}

##########

# sp = '5_Tge'
pdf('figures/anchored_SNPs_and_SVs_overlayed.pdf', width = 10, height = 6)
# pdf(paste0('figures/anchored_SNPs_', sp, '.pdf'), width = 10, height = 6)

for(sp in timemas$codes){
    # tab_filename <- paste0('data/SNP_calls/', sp, '_reduced_filtered_variants.tsv')
    # SV_filename <- paste0('data/SNP_calls/', sp, '_reduced_filtered_variants.tsv')
    # SVs_filt_stringent_union_file <- paste0('data/manta_SV_calls/data/', sp, '/SVs_filt_stringent_union.vcf')
    # mapping_file <- paste0('data/b3v08_anchoring_to_LGs/', sp, '_scf_block_alignment.tsv')
    # mapping_summary_file <- paste0('tables/anchoring/', sp, '_summary.tsv')

    make_variant_density_table_file <- paste0("tables/SNPs/", sp, "_SNPs_on_chromosomes_w", window, ".tsv")

    if ( file.exists(output_file) ){
        variant_density_table <- read.table(make_variant_density_table_file, header = T, sep = '\t', stringsAsFactors = F)
    } else {
        print('please, run D_variant_analysis/make_variant_density_table.R first')
    }

    if ( sp %in% timemas$codes[seq(1, 10, by=2)]){
        plot_SVs_on_LGs(F, T, 'asex', c(0, 1))
    } else {
        plot_SVs_on_LGs(F, T, 'sex', c(0, 22000))
    }

}

dev.off()