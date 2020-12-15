table_name = commandArgs(trailingOnly=TRUE)

source('scripts/R/batch_subset.R')
source('E_assembly_evaluation/batch_table_functions.R')

batches = list()

for(batch in 1:3){
    batches[[batch]] <- batch_stats(batch)
}

batch_table <- data.frame()
for(batch in 1:3){
    batch_table <- rbind(batch_table, get_batch_table(batches[[batch]]))
    write.table(batches[[batch]], paste0('stats/assemblies/batch', batch, '.tsv'), quote = F, sep = '\t', row.names = F)
}

colnames(batch_table) <- c('length_min', 'length_median_asex', 'length_median_sex', 'length_max',
                           'NG50_min', 'NG50_median_asex', 'NG50_median_sex', 'NG50_max',
                           'BUSCOc_min', 'BUSCOc_median_asex', 'BUSCOc_median_sex', 'BUSCOc_max',
                           'BUSCOf_min', 'BUSCOf_median_asex', 'BUSCOf_median_sex', 'BUSCOf_max',
                           'Ns_min', 'Ns_median_asex', 'Ns_median_sex', 'Ns_max')

write.table(batch_table, table_name, quote = F, sep = '\t', row.names = F)
