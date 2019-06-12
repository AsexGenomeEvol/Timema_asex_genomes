library(AsexStats)
source('F_structural_variation_analysis/load_SV_calls.R')
source('F_structural_variation_analysis/filter_SV_calls.R')

# load SV calls
SV_call_manta_files <- paste0("data/", timemas$codes, "/variant_calls/", timemas$codes, "_survivor_manta_calls_union.vcf")
SV_calls <- load_SV_calls(SV_call_manta_files)

get_SV_table_per_ind <- function(heter_sv_table){
    inds = c('00', '01', '02', '03', '04', '05')
    t(sapply(inds, function(ind){
        ind_subset <- heter_sv_table[heter_sv_table[, ind] == 1,]
        table(ind_subset$type)
    }))
}

get_by_type_SV_tables <- function(SV_calls){
    # subset heterozygous and homozygous calls
    heterozygous_SV_calls <- lapply(SV_calls, get_subset_of_SV_calls, 'heterozygous')
    homozygous_SV_calls <- lapply(SV_calls, get_subset_of_SV_calls, 'homozygous')
    polymorphic_SV_calls <- lapply(SV_calls, get_subset_of_SV_calls, 'polymorphic')

    heter_sv_tables <- lapply(heterozygous_SV_calls, get_sv_table)
    homo_sv_tables <- lapply(homozygous_SV_calls, get_sv_table)
    polymorphic_sv_tables <- lapply(polymorphic_SV_calls, get_sv_table)

    # heter_types <- lapply(heter_sv_tables, function(x) { round(table(x$type) / nrow(x), 3) } )
    # heter_types <- lapply(heter_sv_tables, function(x) { table(x$type) } )
    # homo_types <- lapply(homo_sv_tables, function(x) { table(x$type) } )

    heter_sv_summary_tables <- lapply(heter_sv_tables, get_SV_table_per_ind)
    homo_sv_summary_tables <- lapply(homo_sv_tables, get_SV_table_per_ind)
    polymorphic_sv_summary_tables <- lapply(polymorphic_sv_tables, get_SV_table_per_ind)

    merged_sv_summary_tables <- list()
    for( i in 1:10 ){
        heter_tab <- heter_sv_summary_tables[[i]][,1:4]
        homo_tab  <- homo_sv_summary_tables[[i]][,1:4]
        poly_tab   <- polymorphic_sv_summary_tables[[i]][,1:4]

        colnames(heter_tab) <- paste("heter", colnames(heter_tab), sep = "_")
        colnames(homo_tab)  <- paste("homo", colnames(homo_tab), sep = "_")
        colnames(poly_tab)   <- paste("all", colnames(poly_tab), sep = "_")

        correct_order <- as.vector(matrix(1:12, ncol = 4, byrow = T))
        merged_sv_summary_tables[[i]] <- cbind(heter_tab, homo_tab, poly_tab)[,correct_order]
    }

    by_type_tables <- list()
    for ( type in c('DEL', 'DUP', 'INS', 'INV') ){
        # very UGLY soring of all the lists in a single data frame
        type_df <- as.data.frame(matrix(unlist(lapply(merged_sv_summary_tables, function(tab){ tab[, grepl(type, colnames(tab))] } )), nrow = 6))
        if ( ncol(type_df) != 0){
            colnames(type_df) <- paste(rep(timemas$codes, each = 3), rep(c('heter', 'homo', 'all'), 10), sep = "_")
            by_type_tables[[type]] <- type_df
        }
    }

    by_type_tables
}

source('F_structural_variation_analysis/plot_SV_barplots.R')
by_type_tables <- get_by_type_SV_tables(SV_calls)

pdf('figures/SV_overview_zoomed_manta.pdf', width = 10, height = 8)
par(mfrow = c(4,1))
for ( type in c('DEL', 'DUP', 'INS', 'INV') ){
    ymax <- c(5000, 2000, 800, 800)[type == c('DEL', 'DUP', 'INS', 'INV')]
    main <- c('Deletions', 'Duplications', 'Insertions', 'Inversions')[type == c('DEL', 'DUP', 'INS', 'INV')]
    legend <- ifelse(type == 'DEL', T, F)
    plot_barplots(by_type_tables, type, main, ymax, legend)
}
dev.off()

pdf('figures/SV_overview_manta.pdf', width = 20, height = 16)
par(mfrow = c(4,1))
for ( type in c('DEL', 'DUP', 'INS', 'INV') ){
    main <- c('Deletions', 'Duplications', 'Insertions', 'Inversions')[type == c('DEL', 'DUP', 'INS', 'INV')]
    legend <- ifelse(type == 'DEL', T, F)
    plot_barplots(by_type_tables, type, "", NA, legend)
}
dev.off()

SV_call_lumpy_files <- paste0("data/", timemas$codes, "/variant_calls/", timemas$codes, "_survivor_lumpy_calls_union.vcf")
lumpy_SV_calls <- load_SV_calls(SV_call_lumpy_files)
lumpy_by_type_tables <- get_by_type_SV_tables(lumpy_SV_calls)

pdf('figures/SV_overview_lmupy_zoomed.pdf', width = 10, height = 6)
par(mfrow = c(3,1))
for ( type in c('DEL', 'DUP', 'INV') ){
    ymax <- c(5000, 2000, 800)[type == c('DEL', 'DUP', 'INV')]
    main <- c('Deletions', 'Duplications', 'Inversions')[type == c('DEL', 'DUP', 'INV')]
    legend <- ifelse(type == 'DEL', T, F)
    plot_barplots(lumpy_by_type_tables, type, main, ymax, legend)
}
dev.off()

pdf('figures/SV_overview_lmupy.pdf', width = 10, height = 6)
par(mfrow = c(3,1))
for ( type in c('DEL', 'DUP', 'INV') ){
    main <- c('Deletions', 'Duplications', 'Inversions')[type == c('DEL', 'DUP', 'INV')]
    legend <- ifelse(type == 'DEL', T, F)
    plot_barplots(lumpy_by_type_tables, type, main, NA, legend)
}
dev.off()

# Bah, this require a bit of thinking
# How many SVs I expect to find homo/heteroz in all given allelic frequency of the SV?
# get_AA <- function(p){
#     ((1 - p)^2 + 2 * p * (1 - p)) * ((p)^2)^5 +
#     (p^2 + 2 * p * (1 - p)) * ((1 - p)^2)^5
# }
# get_AB <- function(p){
#     2 * p * (1 - p) * ((p)^2)^5 + ...
# }
#
# ps <- seq(0.001, 0.999, length = 1000)
# plot(ps, get_AA(ps), type = 'l', ylim = c(0, 1))
# lines(ps, get_AB(ps), lty = 2)


#### final plot
source('F_structural_variation_analysis/plot_SV_barplots.R')

pdf('figures/SV_overview_manta.pdf', width = 20, height = 16)
par(mfrow = c(4,1))
for ( type in c('DEL', 'DUP', 'INS', 'INV') ){
    plot_barplots(by_type_tables, type, "", NA, F)
}
dev.off()


pdf('figures/SV_merged_overview_manta.pdf', width = 20, height = 4)

source('F_structural_variation_analysis/plot_SV_barplots.R')
merged_types <- list()
merged_types[["MERGED"]] <- by_type_tables[[1]] + by_type_tables[[2]] + by_type_tables[[3]]
plot_barplots(merged_types, "MERGED", "", NA, F)

dev.off()