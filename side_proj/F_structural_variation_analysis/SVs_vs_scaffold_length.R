library(AsexStats)

out_dir <- '~/Desktop/SVs_vs_scf_lengths/'

source('F_structural_variation_analysis/load_SV_calls.R')
source('F_structural_variation_analysis/filter_SV_calls.R')

SV_call_manta_files <- paste0("data/", timemas$codes, "/variant_calls/", timemas$codes, "_survivor_manta_calls_union.vcf")
SV_calls <- load_SV_calls(SV_call_manta_files)
SV_calls <- lapply(SV_calls, get_subset_of_SV_calls, "filter_rare")

tables_of_svs_per_scf <- lapply(SV_calls, function(sp_sv_calls){ table(sapply(sp_sv_calls, function(sv){ sv[1] } )) })

scf_len_files <- paste0("data/", timemas$codes, "/reference/", timemas$codes, "_b3v08_scf.lengths")
scf_lengths <- lapply(scf_len_files, read.table, header = F, col.names = c('scf', 'len'))

SVs_per_scf <- lapply(tables_of_svs_per_scf, function(x){ data.frame(scf = names(x), SVs = as.vector(x)) } )
SVs_per_scf <- lapply(1:10, function(x){ merge(SVs_per_scf[[x]], scf_lengths[[x]]) } )

cum_len_vs_cum_SVs_file <- paste0(out_dir, 'cummulative_SVs_vs_cummulative_lengths.pdf')
pdf(cum_len_vs_cum_SVs_file)
for ( i in 1:10 ){
    plot(   (cumsum(SVs_per_scf[[i]]$len) / 1e6),
            cumsum(SVs_per_scf[[i]]$SVs),
            type = 'l',
            main = timemas$names[[i]],
            xlab = 'cummulative sum of scf lengths [ Mbp ]',
            ylab = 'cumulative sum of number of SVs')
}
dev.off()


svs_vs_cutoff_file <- paste0(out_dir, 'SVs_vs_cutoff.pdf')

pdf(svs_vs_cutoff_file)
for ( i in 1:10 ){
    possible_offs <- as.numeric(names(table(SVs_per_scf[[i]]$len)))
    kept_SVs <- sapply(possible_offs, function(cutoff){ sum(SVs_per_scf[[i]]$SVs[SVs_per_scf[[i]]$len > cutoff]) } )
    plot(   possible_offs,
            kept_SVs,
            type = 'l',
            main = timemas$names[[i]],
            xlim = c(250, 20000),
            xlab = 'cutoff',
            ylab = 'called SVs')
}
dev.off()
