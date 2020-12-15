#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
source('scripts/R/get_optimal_asm.R')
source('E_assembly_evaluation/genome_stat_functions.R')

sp <- args[length(args)]
files <- args[-length(args)]

asm_template <- make_data_frame(variables)
asemblies <- asm_template
for(filename in files){
    asemblies <- rbind(asemblies, make_dataframe_row(filename, asm_template))
}

type <- ifelse(any(grepl('ctg', files)), 'ctgs', 'scfs')
file_to_save <- paste0('stats/assemblies/',sp,'_',type,'.tsv')

if(file.exists(file_to_save)){
    old_asm <- read.table(file_to_save, header = T)
    old_asm$diff_in_sum <- c() # this is only way how to ensure compatibility
    old_asm$score <- c()
    absent_asm <- asemblies[!c(asemblies$dir %in% old_asm$dir),]
    table_to_save <- rbind(old_asm, absent_asm)
    print(paste0('updating: ', file_to_save))
} else {
    print(paste0('creating: ', file_to_save))
    table_to_save <- asemblies
}

# stat widely used
table_to_save$diff_in_sum <- table_to_save$total_sum - 1300000000
table_to_save <- get_asm_scores(table_to_save)
table_to_save <- table_to_save[order(table_to_save$score, decreasing = T),]

write.table(table_to_save, file_to_save, quote=FALSE, sep='\t', row.names = F)
