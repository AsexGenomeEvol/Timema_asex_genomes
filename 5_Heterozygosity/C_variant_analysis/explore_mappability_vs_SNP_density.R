library(AsexStats)
window = 1e6

variant_density_table_files <- paste0("tables/", timemas$codes, "_variants_on_chromosomes_w", window, ".tsv")
variant_density_tables <- lapply(variant_density_table_files, read.table, header = T, sep = '\t', stringsAsFactors = F)

plot(variant_density_tables[[1]]$uniq_mapped ~ variant_density_tables[[1]]$SVs)
