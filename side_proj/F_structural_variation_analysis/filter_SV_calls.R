does_SV_call_carry_a_pattern <- function(x, pattern){
    if ( all(grepl(pattern, x[10:15])) ){
        x
    } else {
        NA
    }
}

is_polymorphic <- function(x, pattern){
    if ( any(grepl("0/1", x[10:15])) && any(grepl("1/1", x[10:15])) ) {
        x
    } else {
        NA
    }
}

is_rare <- function(x, pattern){
    if ( sum(!grepl("./.", x[10:15], fixed = T)) > 1 ) {
        x
    } else {
        NA
    }
}

get_subset_of_SV_calls <- function(sp_sv_calls, sv_subset){
    # the conditional like this allows to set more than one keyword
    if ( 'heterozygous' %in% sv_subset){
        sp_sv_calls <- lapply(sp_sv_calls, does_SV_call_carry_a_pattern, "[.0]/[.1]")
    }
    if ( 'homozygous' %in% sv_subset){
        sp_sv_calls <- lapply(sp_sv_calls, does_SV_call_carry_a_pattern, "[.1]/[.1]")
    }
    if ( 'polymorphic' %in% sv_subset){
        sp_sv_calls <- lapply(sp_sv_calls, is_polymorphic)
    }
    if ( 'filter_rare' %in% sv_subset){
        sp_sv_calls <- lapply(sp_sv_calls, is_rare)
    }
    # space for more sv_subset keywords
    sp_sv_calls[!is.na(sp_sv_calls)]
}

get_SV_table_per_ind <- function(heter_sv_table){
    inds = c('00', '01', '02', '03', '04', '05')
    t(sapply(inds, function(ind){
        ind_subset <- heter_sv_table[heter_sv_table[, ind] == 1,]
        table(ind_subset$type)
    }))
}

get_by_type_SV_tables <- function(SV_calls, filter_rare = F){
    # subset heterozygous and homozygous calls
    if ( filter_rare ){
        het_pattern <- 'heterozygous'
        hom_pattern <- 'homozygous'
        poly_pattern <- 'polymorphic'
    } else {
        het_pattern <- c('heterozygous', 'filter_rare')
        hom_pattern <- c('homozygous', 'filter_rare')
        poly_pattern <- c('polymorphic', 'filter_rare')
    }
    heterozygous_SV_calls <- lapply(SV_calls, get_subset_of_SV_calls, het_pattern)
    homozygous_SV_calls <- lapply(SV_calls, get_subset_of_SV_calls, hom_pattern)
    polymorphic_SV_calls <- lapply(SV_calls, get_subset_of_SV_calls, poly_pattern)

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
