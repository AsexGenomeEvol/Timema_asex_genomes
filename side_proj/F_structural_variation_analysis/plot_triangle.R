plot_trinagle_from_calls <- function(per_type_SV_calls, SV_type, sp){
    i <- which(sp == timemas$codes)
    site_freq <- sapply(per_type_SV_calls[[SV_type]][[i]], alt_alleles)
    pop_heterozygosity <- sapply(per_type_SV_calls[[SV_type]][[i]], gen_heterozygots)
    color = ifelse(i %% 2 == 0, sex_red, asex_blue)

    sfs_file <- paste('figures/sfs/SFS', sp, SV_type,'manta.pdf', sep = '_')

    pdf(sfs_file, width = 9, height = 4)
        barplot(table(site_freq), col = color, cex.axis = 2, cex.names = 1.75, main = paste(sp, SV_type))
    dev.off()

    triangle_file <- paste('figures/triangles/triangle', sp, SV_type,'manta.pdf', sep = '_')

    pdf(triangle_file, width = 6, height = 6)
        plot_trinagle(site_freq, pop_heterozygosity, color, paste(sp, 'heterozygosity vs site frequency'))
    dev.off()
}


plot_trinagle <- function(site_freq, pop_heterozygosity, color, main = 'Triangle'){
    cells <- table(paste(pop_heterozygosity, site_freq))
    total <- sum(cells)
    cells <- round(100 * cells / total, 2)
    cell_coordinates <- lapply(strsplit(names(cells), ' '), as.numeric)
    cr = colorRampPalette(c(rep(color, 4), 'white'))
    k <- kde2d(site_freq, pop_heterozygosity, n = 100, h = 1)

    image(k, col = rev(cr(1000)),
          # ylab = 'number of heterozygous individuals',
          # xlab = 'allele frequency',
          ylim = c(-0.2, 6.2), xlim = c(0.6, 12.4),
          cex.axis = 1.4, bty = 'n', main = main)
    text(sapply(cell_coordinates, function(x){ if( x[2] == 1) { 1.2 } else { x[2] } } ),
         sapply(cell_coordinates, function(x){ x[1] } ),
         cells)
    legend('topright', bty = 'n', paste("Total:", total))
}