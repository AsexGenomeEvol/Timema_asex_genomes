library('RColorBrewer')
library('AsexStats')

pdf('figures/SNP_coverages/one2one_orth_variant_coverages_asex_species.pdf')
# sp = '1_Tps'

get_variant_subset <- function(var_str){
    # var_str <- conserved_pos[,1]
    var_vec <- strsplit(var_str, split = ':')[[1]]
    scf <- var_vec[1]
    range <- as.numeric(strsplit(var_vec[2], split = '-')[[1]])

    variant_tab[variant_tab$scf == scf & variant_tab$pos > range[1] & variant_tab$pos < range[2], ]
}

pal <- brewer.pal(3, 'Dark2')


for (sp in timemas$codes[seq(1,10, by = 2)]){
    tab_filename <- paste0('data/SNP_calls/', sp ,'_reduced_variants.tsv')
    variant_tab <- read.table(tab_filename, stringsAsFactors = F)
    colnames(variant_tab) <- c('scf', 'pos', 'qual', paste0('g', 1:5), paste0('d', 1:5))

    conserved_pos <- read.table(paste0('data/SNP_calls/orth_coords/', substr(sp, 3, 5), '_fivepairs_10sp_orth_coords.txt'), sep = '\t', stringsAsFactors = F)
    table_of_conserved_variants <- do.call(rbind, lapply(conserved_pos[,1], get_variant_subset))

    # hist(variant_tab$d2, breaks = 60)
    # hist(variant_tab$d2[variant_tab$g2 == '1/1'], col = pal[3], add = T, breaks = 60)
    # hist(variant_tab$d2[variant_tab$g2 == '0/0'], col = pal[1], add = T, breaks = 60)
    # hist(variant_tab$d2[variant_tab$g2 == '0/1'], col = pal[2], add = T, breaks = 60)

    ind <- 1
    for(ind in 1:5){
        g <- paste0('g', ind)
        d <- paste0('d', ind)
        name <- paste0(substr(sp, 3, 5), '_', ind)
        # filename <- paste0('figures/SNP_coverages/', name, '_SNP_depth_densities.png')

        if ( any(table_of_conserved_variants[,d] == '.')){
            table_of_conserved_variants[table_of_conserved_variants[,d] == '.', d] <- NA
        }
        table_of_conserved_variants[,d] <- as.numeric(table_of_conserved_variants[,d])

        # hist(variant_tab$d2, breaks = 60, probability = T)
        # png(filename)
            hist(table_of_conserved_variants[table_of_conserved_variants[,g] == '0/0', d], col = pal[1],
                 breaks = 30, probability = T, main = paste(name, 'Distributions of variant depths decomposited by type'),
                 xlab = 'depth')
            hist(table_of_conserved_variants[table_of_conserved_variants[,g] == '1/1', d], col = pal[3], add = T, breaks = 30, probability = T)
            hist(table_of_conserved_variants[table_of_conserved_variants[,g] == '0/1', d], col = pal[2], add = T, breaks = 30, probability = T)
            legend('topright', bty = 'n', pch = 15, col = pal, legend = c('0/0', '0/1', '1/1'))
        # dev.off()
    }
}

dev.off()


depths <- c('d1', 'd2', 'd3', 'd4', 'd5')

for (sp in timemas$codes[seq(1,10, by = 1)]){
    tab_filename <- paste0('data/SNP_calls/', sp ,'_reduced_variants.tsv')
    variant_tab <- read.table(tab_filename, stringsAsFactors = F)
    colnames(variant_tab) <- c('scf', 'pos', 'qual', paste0('g', 1:5), paste0('d', 1:5))

    png(paste0('figures/SNP_qual/', sp ,'_qual_assesment.png'))
        number_of_het_per_variant = rowSums(variant_tab[,c('g1', 'g2', 'g3', 'g4', 'g5')] == '0/1')
        hist(variant_tab$qual[number_of_het_per_variant == 0], breaks = 300, main = sp)
        hist(variant_tab$qual[number_of_het_per_variant == 1], breaks = 300, add = T, col = 'red')
        hist(variant_tab$qual[number_of_het_per_variant > 1], breaks = 300, add = T, col = 'blue')
        lines(c(300, 300), c(0, 1000), col = 'yellow')
    dev.off()

    png(paste0('figures/SNP_qual/', sp ,'_depth_qual_assesment.png'))
        mean_variant_covrage <- rowMeans(apply(variant_tab[, depths], 2, as.numeric))
        plot(mean_variant_covrage[number_of_het_per_variant == 0] ~ variant_tab$qual[number_of_het_per_variant == 0], cex = 0.3, pch = 20, xlab = 'varian quality', ylab = 'mean coveage', main = sp)
        points(mean_variant_covrage[number_of_het_per_variant == 1] ~ variant_tab$qual[number_of_het_per_variant == 1], col = 'red', cex = 0.3, pch = 20)
        points(mean_variant_covrage[number_of_het_per_variant > 1] ~ variant_tab$qual[number_of_het_per_variant > 1], col = 'blue', cex = 0.3, pch = 20)
    dev.off()
}

# hist(variant_tab$qual, breaks = 100)





