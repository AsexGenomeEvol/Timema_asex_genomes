### heterozygosity in intragenomic heterozygosity

# THIS HAVE GONE WRONG SOMEWHERE!!!
# The plots do not make too much sense!!!

library(AsexStats)
source('F_structural_variation_analysis/load_SV_calls.R')
source('F_structural_variation_analysis/plot_SV_barplots.R')
source('F_structural_variation_analysis/filter_SV_calls.R')

# load SV calls
SV_call_manta_files <- paste0("data/", timemas$codes, "/variant_calls/", timemas$codes, "_survivor_manta_calls_union.vcf")
SV_calls <- load_SV_calls(SV_call_manta_files)

sv_genotype_tables <- lapply(SV_calls, get_sv_table, genotypes = T)
# SV_tab <- sv_tables[[1]]

# prepare table
sv_genotype_tables_per_type <- list()
for ( SV_type in c("DEL", "INS", "DUP", "INV")){

    df <- cbind(sapply(sv_genotype_tables, function(x) { apply(x[x$type == SV_type, c(4:9)], 2, table)[2 ,] }),
                sapply(sv_genotype_tables, function(x) { apply(x[x$type == SV_type, c(4:9)], 2, table)[3 ,] }))
    df <- df[, as.vector(matrix(1:20, ncol = 10, byrow = T))]
    colnames(df) <- paste0(rep(timemas$codes, each = 2), c("_het", "_homo"))
    sv_genotype_tables_per_type[[SV_type]] <- df
}

df <- cbind(sapply(sv_genotype_tables, function(x) { apply(x[, c(4:9)], 2, table)[2 ,] }),
            sapply(sv_genotype_tables, function(x) { apply(x[, c(4:9)], 2, table)[3 ,] }))
df <- df[, as.vector(matrix(1:20, ncol = 10, byrow = T))]
colnames(df) <- paste0(rep(timemas$codes, each = 2), c("_het", "_homo"))
sv_genotype_tables_per_type[["ALL"]] <- df

plot_barplots <- function(by_type_tables, SV_type, main = "", ymax = NA, legend = F){
    mins <- apply(by_type_tables[[SV_type]][-1,], 2, min)
    maxes <- apply(by_type_tables[[SV_type]][-1,], 2, max)
    bar_sizes <- apply(by_type_tables[[SV_type]][-1,], 2, median)
    ref_ind <- by_type_tables[[SV_type]][1,]

    ymax <- ifelse(is.na(ymax), max(c(unlist(maxes), unlist(ref_ind))), ymax)
    locations <- barplot(bar_sizes, col = rep(rep(c(asex_blue, sex_red), 10), each = 2),
                         ylim = c(0, ymax), main = main, xaxt = "n", cex.axis=2) # ylab = 'Number of called variants'
    # axis(1, tick=F, labels = , cex.axis=2)
    text(locations[seq(1, by = 2, length = 10)] + 0.5,
         par("usr")[3] - 50, pos = 1,
         xpd = TRUE, labels = timemas$labels, cex = 1.2)

    w <- 0.1
    for ( i in 1:length(mins) ){
        x <- locations[i]
        y_min <- mins[i]
        y_max <- maxes[i]
        lines(c(x, x), c(y_min, y_max), xpd=T, lwd = 1.5)
        lines(c(x - w, x + w), c(y_min, y_min), xpd=T, lwd = 1.5)
        lines(c(x - w, x + w), c(y_max, y_max), xpd=T, lwd = 1.5)
        points(x, ref_ind[i], pch = 20, cex = 1.5)
    }

    if (legend){
        legend('topleft', bty = 'n', pch = 20, 'reference individual')
    }
}

pdf('figures/SV_overview_simple_hetero_homo_manta.pdf', width = 20, height = 16)
    par(mfrow = c(4,1))
    plot_barplots(sv_genotype_tables_per_type, "DEL", "Deletions", legend = T)
    plot_barplots(sv_genotype_tables_per_type, "DUP", "Duplications")
    plot_barplots(sv_genotype_tables_per_type, "INS", "Insertions")
    plot_barplots(sv_genotype_tables_per_type, "INV", "Inversions")
dev.off()

pdf('figures/SV_overview_simple_hetero_homo_all_manta.pdf', width = 20, height = 16)
    plot_barplots(sv_genotype_tables_per_type, "ALL", "", legend = T)
dev.off()
