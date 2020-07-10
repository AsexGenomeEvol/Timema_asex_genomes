library('RColorBrewer')
library('AsexStats')

pdf('figures/SNP_coverages/variant_coverages_all_species.pdf')

for (sp in timemas$codes){
    tab_filename <- paste0('data/SNP_calls/', sp ,'_reduced_variants.tsv')
    variant_tab <- read.table(tab_filename, stringsAsFactors = F)

    colnames(variant_tab) <- c('scf', 'pos', paste0('g', 1:5), paste0('d', 1:5))


    pal <- brewer.pal(3, 'Dark2')

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

        if ( any(variant_tab[,d] == '.')){
            variant_tab[variant_tab[,d] == '.', d] <- NA
            variant_tab[,d] <- as.numeric(variant_tab[,d])
        }

        # hist(variant_tab$d2, breaks = 60, probability = T)
        # png(filename)
            hist(variant_tab[variant_tab[,g] == '0/0', d], col = pal[1],
                 breaks = 60, probability = T, main = paste(name, 'Distributions of variant depths decomposited by type'),
                 xlab = 'depth')
            hist(variant_tab[variant_tab[,g] == '1/1', d], col = pal[3], add = T, breaks = 60, probability = T)
            hist(variant_tab[variant_tab[,g] == '0/1', d], col = pal[2], add = T, breaks = 60, probability = T)
            legend('topright', bty = 'n', pch = 15, col = pal, legend = c('0/0', '0/1', '1/1'))
        # dev.off()
    }
}

dev.off()