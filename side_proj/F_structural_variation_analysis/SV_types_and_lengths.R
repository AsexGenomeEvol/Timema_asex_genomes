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

pdf('figures/SV_sizes.pdf', height = 30, width = 8)
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
    plot_violins_of_one_sv('INV', 'Insertions')
dev.off()

pdf('figures/SV_sizes_duplications_manta.pdf')
    plot_violins_of_one_sv('DUP', 'Duplications')
dev.off()