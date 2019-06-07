library(AsexStats)
source('F_structural_variation_analysis/load_SV_calls.R')
source('F_structural_variation_analysis/filter_SV_calls.R')

# load SV calls
SV_call_manta_files <- paste0("data/", timemas$codes, "/variant_calls/", timemas$codes, "_survivor_manta_calls_union.vcf")
SV_calls <- load_SV_calls(SV_call_manta_files)

# subset heterozygous and homozygous calls
heterozygous_SV_calls <- lapply(SV_calls, get_subset_of_SV_calls, 'heterozygous')
homozygous_SV_calls <- lapply(SV_calls, get_subset_of_SV_calls, 'homozygous')

# load reference info
scf_len_files <- paste0("data/", timemas$codes, "/reference/", timemas$codes, "_b3v08_scf.lengths")
scf_lengths <- lapply(scf_len_files, read.table, header = F, col.names = c('len'), row.names = 1)

# calculate density (occurences / nt, affected nt / nt, statistical ...)
one_sp_calls <- heterozygous_SV_calls[[1]]
one_sp_scfs <- scf_lengths[[1]]

get_density_tables <- function(one_sp_homo_calls, one_sp_hetero_calls, one_sp_scfs){
    one_sp_scfs$homo_SVs <- 0
    one_sp_scfs$homo_SV_lengths <- 0
    one_sp_scfs$hetero_SVs <- 0
    one_sp_scfs$hetero_SV_lengths <- 0

    for (sv in one_sp_homo_calls) {
        one_sp_scfs[sv[1], 'homo_SVs'] = one_sp_scfs[sv[1], 'homo_SVs'] + 1
        sv_len <- as.numeric(ssplit(ssplit(sv[8], split = ";")[3], "=")[2])
        one_sp_scfs[sv[1], 'homo_SV_lengths'] = one_sp_scfs[sv[1], 'homo_SV_lengths'] + sv_len
    }

    for (sv in one_sp_hetero_calls) {
        one_sp_scfs[sv[1], 'hetero_SVs'] = one_sp_scfs[sv[1], 'hetero_SVs'] + 1
        sv_len <- as.numeric(ssplit(ssplit(sv[8], split = ";")[3], "=")[2])
        one_sp_scfs[sv[1], 'hetero_SV_lengths'] = one_sp_scfs[sv[1], 'hetero_SV_lengths'] + sv_len
    }

    one_sp_scfs$hetero_density <- one_sp_scfs$hetero_SVs / one_sp_scfs$len
    one_sp_scfs$homo_density <- one_sp_scfs$homo_SVs / one_sp_scfs$len
    one_sp_scfs
}

scf_densities <- list()
for (i in 1:10) {
    scf_densities[[i]] <- get_density_tables(homozygous_SV_calls[[i]], heterozygous_SV_calls[[i]], scf_lengths[[i]])
}

#### This would be better in python



#
# mean(hetero_density == 0)
# mean(homo_density == 0)
#
# non_zero_scfs <- hetero_density != 0 | homo_density != 0
#
# plot(hetero_density[non_zero_scfs] ~ homo_density[non_zero_scfs], pch = 20)