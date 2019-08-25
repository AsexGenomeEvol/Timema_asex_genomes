library(AsexStats)
source('F_structural_variation_analysis/load_SV_calls.R')
source('F_structural_variation_analysis/filter_SV_calls.R')

# load SV calls
SV_call_manta_files <- paste0("data/", timemas$codes, "/variant_calls/", timemas$codes, "_survivor_manta_calls_union.vcf")
SV_calls <- load_SV_calls(SV_call_manta_files)

# subset heterozygous and homozygous calls
heterozygous_SV_calls <- lapply(SV_calls, get_subset_of_SV_calls, 'heterozygous')
homozygous_SV_calls <- lapply(SV_calls, get_subset_of_SV_calls, 'homozygous')

heter_sv_tables <- lapply(heterozygous_SV_calls, get_sv_table)
homo_sv_tables <- lapply(homozygous_SV_calls, get_sv_table)
all_sv_tables <- lapply(SV_calls, get_sv_table)

# heter_types <- lapply(heter_sv_tables, function(x) { round(table(x$type) / nrow(x), 3) } )
heter_types <- lapply(heter_sv_tables, function(x) { table(x$type) } )
homo_types <- lapply(homo_sv_tables, function(x) { table(x$type) } )

# something like this to per SV plotting
# par(mfrow = c(2,5))
# ylim <- c(0, max(unlist(heter_types)))
# for(i in c(1,3,5,7,9)){
#     barplot(heter_types[[i]], col = asex_blue, main = timemas$labels[i], ylim = ylim, cex.axis = 1.6, cex.names = 1.6)
# }
# for(i in c(2,4,6,8,10)){
#     barplot(heter_types[[i]], col = sex_red, main = timemas$labels[i], ylim = ylim, cex.axis = 1.6, cex.names = 1.6)
# }

# length plotting

require(digest)
library(sm)
source('F_structural_variation_analysis/plot_violins.R')

pdf('figures/SV_sizes_manta.pdf', height = 30, width = 8)
par(mfrow = c(5,1))
for (i in seq(1, 10, by = 2)){
    label <- timema_pairs$labels[(i + 1) / 2]
    asex_sp <- all_sv_tables[[i]]
    asex_sp$sex <- 'asex'
    sex_sp <- all_sv_tables[[i + 1]]
    sex_sp$sex <- 'sex'
    one_pair <- rbind(asex_sp, sex_sp)
    one_pair$loglen <- log10(one_pair$len)
    boxplot(loglen ~ sex + type, one_pair, col = c(asex_blue, sex_red), main = label)
}
dev.off()

plot_violins_of_one_sv <- function(type, main){
    sv_sizes <- lapply(all_sv_tables, function(x) { x$len[x$type == type] } )
    plotLogViolins(sv_sizes)
    axis(2, at = c(1, 2, 3, 4, 5), labels = c('10', '100', '1000', '10000', '100000'), cex.axis = 1)
    title(main)
}

# DEL   DUP   INS   INV   TRA
pdf('figures/SV_sizes_deletions_manta.pdf')
    plot_violins_of_one_sv('DEL', 'Deletions')
dev.off()

pdf('figures/SV_sizes_inversions_manta.pdf')
    plot_violins_of_one_sv('INV', 'Inversions')
dev.off()

pdf('figures/SV_sizes_insertions_manta.pdf')
    plot_violins_of_one_sv('INS', 'Insertions')
dev.off()

pdf('figures/SV_sizes_duplications_manta.pdf')
    plot_violins_of_one_sv('DUP', 'Duplications')
dev.off()


##### PLOT IT PER HOMO/HETERO

# get non-rare
homo_sv_tables <- lapply(homo_sv_tables,
                         function(one_sp){ one_sp[apply(one_sp[,4:9], 1, sum ) > 1 & one_sp[,4] == 0, ]} )

heter_sv_tables <- lapply(heter_sv_tables,
                         function(one_sp){ one_sp[apply(one_sp[,4:9], 1, sum ) > 1, ]} )

# one_sp <- homo_sv_tables[[1]]
for (type_to_plot in c('INV', 'INS', 'DEL')){
    pdf(paste0('figures/SV_hetero_homo_sizes_', type_to_plot , '_manta.pdf'), width = 8, height = 8)

    homo_sizes <- lapply(homo_sv_tables, function(x){ x[x$type == type_to_plot, 'len'] } )
    hetero_sizes <- lapply(heter_sv_tables, function(x){ x[x$type == type_to_plot, 'len'] } )

    ylim <- c(min(log10(abs(unlist(c(homo_sizes, hetero_sizes))))),
              max(log10(abs(unlist(c(homo_sizes, hetero_sizes))))))
    par(mfrow = c(2, 1))
    # ASEX
    asexuals <- seq(1, 10, by = 2)
    par(mar = c(3.2, 3.5, 2.5, 1))
    plot(x=NULL, y=NULL,
         xlim = c(0.5, 5.5), ylim = ylim,
         type="n", ann=FALSE, axes=F)
    axis(1, at=c(1:5), labels = F)
    axis(2, at = c(1, 2, 3, 4, 5), labels = c('10', '100', '1000', '10000', '100000'), cex.axis = 1)
    text(c(1:5), par("usr")[3] - 0.3, srt = 15, pos = 1, xpd = TRUE, labels = timemas$labels[asexuals], cex = 1.3)
    mtext(expression(paste("Sizes [" , log[10], " nt]")), side = 2, line = +2, cex = 1.6)
    title(type_to_plot)
    legend('topright', c('left - all homozygous', 'right - all heterozygous'), bty = 'n')

    for(i in 1:5){
        sp_to_plot <- asexuals[i]
        vioplot2(log10(abs(homo_sizes[[sp_to_plot]])),
                 at = i, side = "left", col = asex_blue,
                 add = T, h = 0.2)
        vioplot2(log10(abs(hetero_sizes[[sp_to_plot]])),
                 at = i, side = "right", col = asex_blue,
                 add = T, h = 0.2)
    }

    par(mar = c(3.2, 3.5, 1, 1))
    sexuals <- seq(2, 10, by = 2)
    plot(x=NULL, y=NULL,
         xlim = c(0.5, 5.5), ylim = ylim,
         type="n", ann=FALSE, axes=F)
    axis(1, at=c(1:5), labels = F)
    axis(2, at = c(1, 2, 3, 4, 5), labels = c('10', '100', '1000', '10000', '100000'), cex.axis = 1)
    text(c(1:5), par("usr")[3] - 0.3, srt = 15, pos = 1, xpd = TRUE, labels = timemas$labels[sexuals], cex = 1.3)
    mtext(expression(paste("Sizes [" , log[10], " nt]")), side = 2, line = +2, cex = 1.6)

    for(i in 1:5){
        sp_to_plot <- sexuals[i]
        vioplot2(log10(abs(homo_sizes[[sp_to_plot]])),
                 at = i, side = "left", col = sex_red,
                 add = T, h = 0.2)
        vioplot2(log10(abs(hetero_sizes[[sp_to_plot]])),
                 at = i, side = "right", col = sex_red,
                 add = T, h = 0.2)
    }
    dev.off()
}
