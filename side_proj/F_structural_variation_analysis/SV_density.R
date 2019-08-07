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
scf_lengths <- lapply(scf_len_files, read.table, header = F, col.names = c('scf', 'len'), row.names = 1)

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

# pause
# for (i in 1:10) {
#     write.table(scf_densities[[i]], paste0("stats/", timemas$codes[i], "_SV_densities_manta.tsv"), quote = F, sep = '\t')
# }

# unpause
# scf_densities <- list()
# for (i in 1:10) {
#     scf_densities[[i]] <- read.table(paste0("stats/", timemas$codes[i], "_SV_densities_manta.tsv"))
# }

get_expectation_tables <- function(one_sp){
    one_sp <- one_sp[one_sp$len > 100000,]

    overall_homo_sv_den <- sum(one_sp$homo_SVs) / sum(one_sp$len)
    overall_hetero_sv_den <- sum(one_sp$hetero_SVs) / sum(one_sp$len)
    one_sp$homo_SV_probs <- apply(one_sp, 1, function(x){ pbinom(x[2], x[1], overall_homo_sv_den) } )
    one_sp$hetero_SV_probs <- apply(one_sp, 1, function(x){ pbinom(x[4], x[1], overall_hetero_sv_den) } )
    one_sp$expected_homo_SVs <- apply(one_sp, 1, function(x){ x[1] * overall_homo_sv_den } )
    one_sp$expected_hetero_SVs <- apply(one_sp, 1, function(x){ x[1] * overall_hetero_sv_den } )
    one_sp$hetero_SVs_dif <- one_sp$hetero_SVs - one_sp$expected_hetero_SVs
    one_sp$homo_SVs_dif <- one_sp$homo_SVs - one_sp$expected_homo_SVs
    one_sp
}

SV_expecations <- lapply(scf_densities, get_expectation_tables)

## Does not lead anywhere. Not sure how comes that the last time I had the |_ plot. (shared by 2?)
## Need to rething it
## Possible courses of actions are
### - mapping our scaffolds to the reference asm
### - filtering of low quality SVs to reduce noise
# lapply(SV_expecations, function(x) { cor.test(x$hetero_density, x$homo_density) } )
# lapply(SV_expecations, function(x) { cor.test(x$hetero_SVs, x$homo_SVs) } )
# # # pbinom(one_sp$homo_SVs[1], one_sp$len[1], overall_homo_sv_den)
# #
# # one_sp <- one_sp[one_sp$expected_homo_SVs > 0.5,]
# boxplot(one_sp$expected_homo_SVs ~ one_sp$homo_SVs, xlab = "scaffolds given observed SVs", ylab = "distributions of expectations of scaffolds")
# lines(c(0,12), c(0, 12))
#
# hist(one_sp$hetero_SVs_dif, col = 'purple')
# hist(one_sp$homo_SVs_dif, col = 'orange')
#
# plot(one_sp$hetero_SVs_dif ~ one_sp$homo_SVs_dif)
#
# library(MASS)
# library(RColorBrewer)
# rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
#
# k <- kde2d(one_sp$hetero_SVs, one_sp$homo_SVs, n=100)
# k$z <- log10(k$z)
# image(k, col=rf(20)[1:15])
#       xlab = 'number of heterozygous individuals',
#       ylab = 'allele frequency')
