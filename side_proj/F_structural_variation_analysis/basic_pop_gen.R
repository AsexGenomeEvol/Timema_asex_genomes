library(AsexStats)
source('F_structural_variation_analysis/plot_SFS.R')
source('F_structural_variation_analysis/load_SV_calls.R')

###############
# MANTA CALLS #
###############

SV_call_manta_files <- paste0("data/", timemas$codes, "/variant_calls/", timemas$codes, "_survivor_manta_calls_union.vcf")
SV_calls <- load_SV_calls(SV_call_manta_files)

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

SV_call_lumpy_files <- paste0("data/", timemas$codes, "/variant_calls/", timemas$codes, "_survivor_lumpy_calls_union.vcf")
SV_calls <- load_SV_calls(SV_call_lumpy_files)
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
