sample_table <- read.table('tables/resequencing_samples.tsv', header = T)

get_coverage <- function(x){
    sp <- x['sp_ID']
    sample <- x['sample_ID']
    cov_file <- paste0('data/', sp, '/trimmed_reads/', sample, '/', sp, '_coverage.stats')
    coverages <- read.table(cov_file)$V2
    return( sum(coverages) / 1.381e9 )
}

sample_table$cov <- round(apply(sample_table, 1, get_coverage), 1)

write.table(sample_table, 'A_mapping_and_preprocessing/resequencing_samples',
            sep = '\t', quote = F, row.names = F)
