library(AsexStats)
library(RColorBrewer)

###############################
##### plotting parameters #####
###############################

window = 1e6
gap_beween_chromosomes = 3

###############################

get_chr_size <- function(chr){
    # when scaffolds are mapped to reference, I always leave 10000 bases between them on a linkage group (consistent with the NCBI reference that used this arbitrary number of bases)
    sum(reference[reference$chromosome == chr,'len']) + ((sum(reference$chromosome == chr) - 1) * 10000)
}

plot_SVs_on_LGs <- function(exclude_zero = T, lines = F, pal = 'asex'){

    # BrBG
    if ( pal == 'asex' ){
        pal <- brewer.pal(5, "YlGnBu")[c(3,5)]
    } else {
        pal <- brewer.pal(5, "YlOrRd")[c(3,5)]
    }
    # pal <- addalpha(pal, alpha)

    ylim = range(variant_density_table$variants, na.rm = T)
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
        variant_density_table[variant_density_table$common == 0, 'common'] <- NA
    }

    if ( lines ){
        lines(variant_density_table$variants, col = pal[1], lwd = 1.6)
        lines(variant_density_table$common, col = pal[2], lwd = 1.6)
    } else {
        points(variant_density_table$variants, pch = 20, col = pal[1])
        points(variant_density_table$common, pch = 20, col = pal[2])
    }
    legend('topright', bty = 'n', c('all variants', 'common'), pch = 20, col = pal[c(1,2)])
}

get_lg_windows <- function(i) {
    all_windows <- seq(0, chromosomes[i,'rounded_len'], by = window)
    adj <- chromosomes[i, 'adjustments']
    lg_from <- all_windows[1:(length(all_windows) - 1)]
    lg_to <- all_windows[2:length(all_windows)]
    data.frame(lg = chromosomes[i, 'chr'],
               lg_from = lg_from,
               lg_to = lg_to,
               genome_from = lg_from + adj,
               genome_to = lg_from + adj)
}

### LOADING INFO ABOUT GENOMES
# scf_len_file <- paste0("data/", sp, "/reference/", sp, "_b3v08_scf.lengths")
# scf_lengths <- read.table(scf_len_file, header = F, row.names = 1, col.names = c('scf', 'len'))
# scf_to_analyze <- lapply(scf_lengths, function(x){rownames(x)[x[,1] > 100000]})
# > head(scf_lengths)
#                            len
# 5_Tge_b3v08_scaf000001 1405450
# 5_Tge_b3v08_scaf000002 1258996
# source('../timema_assembly/scripts/prepare_reference.R')

reference <- read.table('data/external_ref/sex_lg_assigment_scores_1.4a.tsv')
colnames(reference) <- c('scf_o', 'scf', 'score', 'cov', 'len', 'asignment')
reference$chromosome <- sapply(strsplit(reference$scf, "_"), function(x) { x[1] } )

chromosomes <- data.frame(chr = paste0('lg', c(1:12, 'X')))
chromosomes$len <- sapply( chromosomes$chr, get_chr_size )

##########

# sp = '5_Tge'

### LOADING SNPs
pdf('figures/anchored_SNPs.pdf', width = 10, height = 6)
# pdf(paste0('figures/anchored_SNPs_', sp, '.pdf'), width = 10, height = 6)

for(sp in timemas$codes){
    tab_filename <- paste0('data/SNP_calls/', sp, '_reduced_filtered_variants.tsv')
    variant_tab <- read.table(tab_filename, stringsAsFactors = F)
    colnames(variant_tab) <- c('scf', 'pos', 'qual', paste0('g', 1:5), paste0('d', 1:5), 'ref_scf', 'ref_pos', 'lg', 'lg_pos')

    if ( sp == '1_Tps'){
        variant_tab[,'g4'] <- NA
        variant_tab[,'g5'] <- NA
    }
    if ( sp == '4_Tte'){
        variant_tab[,'g1'] <- NA
    }
    if ( sp == '2_Tsi'){
        variant_tab[,'g1'] <- NA
        variant_tab[,'g2'] <- NA
        variant_tab[,'g3'] <- NA
    }

    anchored_variant_tab <- variant_tab[!is.na(variant_tab$lg_pos), ]

    # I need to round up chromosome lengths to avoid one window being part of two chromosomes
    chromosomes$rounded_len <- ceiling(chromosomes$len / window) * window + (window * gap_beween_chromosomes)
    chromosomes$adjustments <- cumsum(c(0, chromosomes$rounded_len[1:(nrow(chromosomes) - 1)]))
    rownames(chromosomes) <- chromosomes$chr
    chromosomes <- chromosomes[chromosomes$chr != 'lgX', ]

    #### prepare variant template
    variant_density_table <- do.call("rbind", lapply(1:12, get_lg_windows))
    variant_density_table$variants <- 0
    variant_density_table$common <- 0

    for (i in 1:nrow(anchored_variant_tab)) {
        lg = anchored_variant_tab[i, 'lg']
        pos = anchored_variant_tab[i, 'lg_pos']
        row <- variant_density_table$lg == lg & variant_density_table$lg_from < pos & variant_density_table$lg_to > pos
        variant_density_table[row, 'variants'] <- variant_density_table[row, 'variants'] + 1
        derived_variants <- sum(anchored_variant_tab[i, paste0('g', 1:5)] == '1/1', na.rm = T)
        derived_variants_freq <- mean(anchored_variant_tab[i, paste0('g', 1:5)] == '1/1', na.rm = T)
        if ( derived_variants == 0){
            next
        }
        if ( derived_variants_freq != 1 & derived_variants != 1 ){
            variant_density_table[row, 'common'] <- variant_density_table[row, 'common'] + 1
        }
    }

    if ( sp %in% timemas$codes[seq(1, 10, by=2)]){
        plot_SVs_on_LGs(F, T, 'asex')
    } else {
        plot_SVs_on_LGs(F, T, 'sex')
    }

}

dev.off()