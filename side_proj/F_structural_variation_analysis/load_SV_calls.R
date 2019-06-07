# load a species list of SVs lists

load_SV_calls <- function(files, headers = F){
    SV_calls <- lapply(files, readLines)
    SV_calls <- lapply(SV_calls, function(x) { x[!grepl("^##", x)]} )
    if ( headers ){
        lapply(SV_calls, function(x) { unlist(strsplit(x[1], '\t'))} )
    } else {
        lapply(SV_calls, function(x) { strsplit(x[-1], '\t')} )
    }
}

get_sv_table <- function(one_sp_sv_calls){
     info_vec <- sapply(one_sp_sv_calls, function(x){x[8]} )
     info_vec <- strsplit(info_vec, ';')
     data.frame(
          scf = sapply(info_vec, get_entry, "CHR2"),
          len = as.numeric(sapply(info_vec, get_entry, "AVGLEN")),
          type = sapply(info_vec, get_entry, "SVTYPE"),
          presence = sapply(info_vec, get_entry, "SUPP_VEC")
     )
}

get_entry <- function(infoline, entry_name){
    entry <- ssplit(grep(entry_name, infoline, value = T))
    ifelse(length(entry) == 0, NA, entry[2])
}