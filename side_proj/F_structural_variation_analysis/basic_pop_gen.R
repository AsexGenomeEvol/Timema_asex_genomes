library(AsexStats)
source('F_structural_variation_analysis/plot_SFS.R')

###############
# MANTA CALLS #
###############

SV_call_files <- paste0("data/", timemas$codes, "/variant_calls/", timemas$codes, "_survivor_manta_calls_union.vcf")

SV_calls <- lapply(SV_call_files, readLines)
SV_calls <- lapply(SV_calls, function(x) { x[!grepl("^##", x)]} )
SV_headers <- lapply(SV_calls, function(x) { unlist(strsplit(x[1], '\t'))} )
SV_calls <- lapply(SV_calls, function(x) { strsplit(x[-1], '\t')} )
# species list of SVs lists

# SFSs
pdf('figures/SFS_manta.pdf', width = 12, height = 8)
layout(mat = matrix(1:10, nrow = 2, ncol = 5))
for( i in 1:10){
    plot_SFS(SV_calls[[i]], timemas$labels[i], timemas$cols[i])
}
dev.off()

pdf('figures/SFS_vs_heterozygosity_manta.pdf', width = 30, height = 12)
layout(mat = matrix(1:10, nrow = 2, ncol = 5))
for( i in 1:10){
    plot_SFS_vs_heterozygosity(SV_calls[[i]], timemas$labels[i])
}
dev.off()


################
# SMOOVE CALLS #
################


SV_call_files <- paste0("data/", timemas$codes, "/variant_calls/", timemas$codes, "_survivor_lumpy_calls_union.vcf")

SV_calls <- lapply(SV_call_files, readLines)
SV_calls <- lapply(SV_calls, function(x) { x[!grepl("^##", x)]} )
SV_headers <- lapply(SV_calls, function(x) { unlist(strsplit(x[1], '\t'))} )
SV_calls <- lapply(SV_calls, function(x) { strsplit(x[-1], '\t')} )
# species list of SVs lists

# SFSs
pdf('figures/SFS_smove.pdf', width = 12, height = 8)
layout(mat = matrix(1:10, nrow = 2, ncol = 5))
for( i in 1:10){
    plot_SFS(SV_calls[[i]], timemas$labels[i], timemas$cols[i])
}
dev.off()

pdf('figures/SFS_vs_heterozygosity_smove.pdf', width = 30, height = 12)
layout(mat = matrix(1:10, nrow = 2, ncol = 5))
for( i in 1:10){
    plot_SFS_vs_heterozygosity(SV_calls[[i]], timemas$labels[i])
}
dev.off()
