library(AsexStats)
source('F_structural_variation_analysis/plot_SFS.R')
source('F_structural_variation_analysis/load_SV_calls.R')
source('F_structural_variation_analysis/plot_triangle.R')

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

##################
# PLOT TRIANGLES #
##################

for ( SV_type in c("DEL", "INS", "DUP", "INV")){
    for ( sp in timemas$codes ) {
        plot_trinagle_from_calls(per_type_SV_calls, SV_type, sp)
    }
}

