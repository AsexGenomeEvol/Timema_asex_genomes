library(AsexStats)

source('B2_genome_profiling/extract_heterozygosity_genomescope.R')
genomescope_files <- sapply(timemas$codes, function(x){(dir(paste0('data/',x,'/genomescope'), pattern = 'summary', full.names = T))})

heterozygosities_min <- sapply(genomescope_files, extract_heterozygosity, min = T, rounding = 6)
heterozygosities_max <- sapply(genomescope_files, extract_heterozygosity, min = F, rounding = 6)
heterozygosities <- (heterozygosities_min + heterozygosities_max) / 2

genome_tab$heterozygosity <- heterozygosities

png("figures/heterozygosity_continuity_relation.png")

plot(genome_tab$N50 ~ genome_tab$heterozygosity, pch = 20, col = c(asex_blue, sex_red),
     ylab = 'Wigthed median of scaffold sizes N50 [ kbp ]', xlab = ' Heterozygosity [ % ] ', main = 'Lower continuity of assembly explained by elevated heterozygosity')

dev.off()