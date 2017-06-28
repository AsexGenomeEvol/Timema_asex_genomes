getSVTYPE <- function(info){
    ssplit(grep("SVTYPE", info, value = T))[2]
}

getSVLEN <- function(info){
    SVLEN <- ssplit(grep("SVLEN", info, value = T))
    ifelse(length(SVLEN) == 0, NA, as.numeric(SVLEN[2]))
}

load_sv <- function(sp){
    file <- paste0('data/', sp , '/variant_calls/ref_is350/manta/results/variants/diploidSV.vcf')
    vcf_sv_table <- read.table(file, stringsAsFactors = F)
    colnames(vcf_sv_table) <- SV_colnames
    # vcf_sv_table <- vcf_sv_table[vcf_sv_table$FILTER == 'PASS',]
    vcf_sv_table_info <- strsplit(vcf_sv_table$INFO, ';')
    vcf_sv_table$SVTYPE <- unlist(lapply(vcf_sv_table_info, getSVTYPE))
    vcf_sv_table$SVLEN <- unlist(lapply(vcf_sv_table_info, getSVLEN))
    return(vcf_sv_table)
}
