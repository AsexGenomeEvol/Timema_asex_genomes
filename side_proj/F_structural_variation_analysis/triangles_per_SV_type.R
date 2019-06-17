library(AsexStats)
source('F_structural_variation_analysis/plot_SFS.R')
source('F_structural_variation_analysis/load_SV_calls.R')

###############
# MANTA CALLS #
###############

SV_call_manta_files <- paste0("data/", timemas$codes, "/variant_calls/", timemas$codes, "_survivor_manta_calls_union.vcf")
SV_calls <- load_SV_calls(SV_call_manta_files)

list_of_types <- lapply( SV_calls, function(x){
    info_vec <- sapply(x, function(sv_call) { sv_call[8] })
    info_vec <- strsplit(info_vec, ';')
    sapply(info_vec, get_entry, "SVTYPE")} )

per_type_SV_calls <- list()
for ( SV_type in c("DEL", "INS", "DUP", "INV")){
    list_of_subsests <- lapply(list_of_types, function(x) { which(x == SV_type) } )
    filtered_SV_calls <- lapply(1:10, function(sp_index){ SV_calls[[sp_index]][list_of_subsests[[sp_index]]] } )
    per_type_SV_calls[[SV_type]] <- filtered_SV_calls
}

plot_trinagle <- function(per_type_SV_calls, SV_type, sp){
    i <- which(sp == timemas$codes)
    site_freq <- sapply(per_type_SV_calls[[SV_type]][[i]], alt_alleles)
    pop_heterozygosity <- sapply(per_type_SV_calls[[SV_type]][[i]], gen_heterozygots)

    sfs_file <- paste('figures/sfs/SFS', sp, SV_type,'manta.pdf', sep = '_')

    pdf(sfs_file, width = 9, height = 4)
        barplot(table(site_freq), col = ifelse(i %% 2 == 0, sex_red, asex_blue), cex.axis = 2, cex.names = 1.75, main = paste(sp, SV_type))
    dev.off()

    cells <- table(paste(pop_heterozygosity, site_freq))
    total <- sum(cells)
    cells <- round(100 * cells / total, 2)
    cell_coordinates <- lapply(strsplit(names(cells), ' '), as.numeric)
    if ( i %% 2 == 0 ){
        cr = colorRampPalette(c(sex_red, sex_red, sex_red, sex_red, 'white'))
    } else {
        cr = colorRampPalette(c(asex_blue, asex_blue, asex_blue, asex_blue, 'white'))
    }
    k <- kde2d(Tce_site_freq, Tce_pop_heterozygosity, n = 100)

    triangle_file <- paste('figures/triangles/triangle', sp, SV_type,'manta.pdf', sep = '_')

    pdf(triangle_file, width = 6, height = 6)
        image(k, col = rev(cr(1000)),
              # ylab = 'number of heterozygous individuals',
              # xlab = 'allele frequency',
              ylim = c(-0.2, 6.2), xlim = c(0.6, 12.4),
              cex.axis = 1.4, bty = 'n', main = paste(sp, SV_type))
        text(sapply(cell_coordinates, function(x){ if( x[2] == 1) { 1.2 } else { x[2] } } ),
             sapply(cell_coordinates, function(x){ x[1] } ),
             cells)
        legend('topright', bty = 'n', paste("Total:", total))
    dev.off()
}

for ( SV_type in c("DEL", "INS", "DUP", "INV")){
    for ( sp in timemas$codes ) {
        plot_trinagle(per_type_SV_calls, SV_type, sp)
    }
}

