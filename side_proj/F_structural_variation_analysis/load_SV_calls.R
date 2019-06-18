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

get_sv_table <- function(one_sp_sv_calls, genotypes = F, coordinates = F){
     info_vec <- sapply(one_sp_sv_calls, function(x){x[8]} )
     info_vec <- strsplit(info_vec, ';')
     sv_frame <- data.frame(
          scf = sapply(info_vec, get_entry, "CHR2"),
          len = as.numeric(sapply(info_vec, get_entry, "AVGLEN")),
          type = sapply(info_vec, get_entry, "SVTYPE")
     )
     if (genotypes) {
         presence_vec = sapply(one_sp_sv_calls, get_genotypes)
         sv_frame[, c('00', '01', '02', '03', '04', '05')] <- matrix(presence_vec, ncol = 6, byrow = T)
     } else {
         presence_vec = as.numeric(ssplit(sapply(info_vec, get_entry, "SUPP_VEC"), split = ''))
         sv_frame[, c('00', '01', '02', '03', '04', '05')] <- matrix(presence_vec, ncol = 6, byrow = T)
     }
     if (coordinates) {
         sv_frame$pos <- as.numeric(sapply(one_sp_sv_calls, function(x) { x[2] } ))
         sv_frame$end <- as.numeric(sapply(info_vec, get_entry, "END"))
     }
     sv_frame
}

get_entry <- function(infoline, entry_name){
    entry <- ssplit(grep(entry_name, infoline, value = T))
    ifelse(length(entry) == 0, NA, entry[2])
}

get_genotypes <- function(infoline){
    substr(infoline[10:15], 1, 3)
}