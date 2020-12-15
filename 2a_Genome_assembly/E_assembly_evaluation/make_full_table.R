asm_type = commandArgs(trailingOnly=TRUE)

pattern <- paste0('*',asm_type,'.tsv') # '*ctgs.tsv'
asm_stats <- data.frame()
contig_files <- dir('assemblies', pattern = pattern, full.names = T)
species <- substr(dir('assemblies', pattern = pattern), 1, 5)

for(i in 1:length(species)){
   tmp_tab <- read.table(contig_files[i], header = T)
   tmp_tab$sp <- species[i]
   tmp_tab$diff_in_sum <- c()
   asm_stats <- rbind(asm_stats, tmp_tab)
}
asm_stats$soft <- factor(asm_stats$soft, levels = c("SOAP", "abyss", "BESST"))
asm_stats$diff_in_sum <- asm_stats$total_sum - 1300000000


write.table(asm_stats,
            paste0('assemblies/', asm_type, '_fulltable.tsv'),
            quote = F, sep = '\t', row.names = F)
