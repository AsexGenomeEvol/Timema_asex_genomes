does_SV_call_carry_a_pattern <- function(x, pattern){
    if ( all(grepl(pattern, x[10:15])) ){
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
    # space for more sv_subset keywords
    sp_sv_calls[!is.na(sp_sv_calls)]
}