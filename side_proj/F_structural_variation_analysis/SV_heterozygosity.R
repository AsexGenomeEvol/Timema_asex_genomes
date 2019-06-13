library(AsexStats)
source('F_structural_variation_analysis/load_SV_calls.R')
source('F_structural_variation_analysis/filter_SV_calls.R')
source('F_structural_variation_analysis/plot_SV_barplots.R')

SV_call_manta_files <- paste0("data/", timemas$codes, "/variant_calls/", timemas$codes, "_survivor_manta_calls_union.vcf")
SV_calls <- load_SV_calls(SV_call_manta_files)

get_number_of_heterozygous_loci_in_reference <- function(one_sp_sv_calls){
     info_vec <- sapply(one_sp_sv_calls, function(x){x[10]} )
     info_vec <- substr(info_vec, 1, 3)
     sum(info_vec == '0/1')
}

heterozygous_SVs <- sapply(SV_calls, get_number_of_heterozygous_loci_in_reference)

# To consider -> use all 6 individuals for heterozygosity estimate, not just one
# that would allow ~ CI

pdf('figures/SV_heterozygosity.pdf')
     locations <- barplot(heterozygous_SVs, col = c(asex_blue, sex_red), cex.axis = 1.8)
     text(locations - 0.3, par("usr")[3] - ((max(heterozygous_SVs) - min(heterozygous_SVs)) / 40), srt = 20, pos = 1, xpd = TRUE, labels = timemas$labels, cex = 1.2)
# , col = c(asex_blue, sex_red)
dev.off()