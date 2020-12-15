require(AsexStats)

read_busco <- function(asm_path){
    if(is.na(asm_path)){
        return(data.frame(complete = NA, duplicated = NA, fragmented = NA, missing = NA))
    }
    busco_file <- dir(asm_path, full.names = T, pattern = 'short_summary')
    if(length(busco_file) != 1){
        return(data.frame(complete = NA, duplicated = NA, fragmented = NA, missing = NA))
    }
    busco_file <- readLines(busco_file)
    total_genes <- as.numeric(ssplit(busco_file[15], '\t')[2])
    bscores <- data.frame(complete = as.numeric(ssplit(busco_file[10], '\t')[2]),
                          duplicated = as.numeric(ssplit(busco_file[12], '\t')[2]),
                          fragmented = as.numeric(ssplit(busco_file[13], '\t')[2]),
                          missing = as.numeric(ssplit(busco_file[14], '\t')[2]))
    bscores[1,] <- round(100 * (bscores[1,] / total_genes), 2)
    return(bscores)
}

batch_stats <- function(batch){
    batch_table <- make_data_frame(c('sp', 'total_sum', 'NG50', 'BUSCOc', 'BUSCOf', 'Ns'))
    for(sp in timemas){
        asm <- read.table(paste0('stats/assemblies/',sp,'_scfs.tsv'), header = T)
        sp_batch <- batch_subset(asm, batch)[1,]
        asm_path <- paste0('data/',sp,'/assembly/',sp_batch$dir)
        N_file <- dir(asm_path, full.names = T, pattern = '_Ns.tsv')
        if(length(N_file) == 1){
            Ns <- read.table(N_file)[1,2]
        } else {
            Ns <- NA
        }
        BUSCOs <- read_busco(asm_path)
        batch_table <- rbind(batch_table, data.frame(sp = sp,
                                                     total_sum = sp_batch$total_sum,
                                                     NG50 = sp_batch$NG50,
                                                     complete = BUSCOs$complete,
                                                     fragmented = BUSCOs$fragmented,
                                                     Ns = Ns))
    }
    batch_table$Ns <- round(100 * (batch_table$Ns / batch_table$total_sum), 2)
    batch_table$total_sum <- round(batch_table$total_sum / 1e9, 3)
    return(batch_table)
}

var_summary <- function(batch_table, var){
    return(c(min(batch_table[,var], na.rm = T),
             median(batch_table[seq(1,10, by = 2),var], na.rm = T),
             median(batch_table[seq(2,10, by = 2),var], na.rm = T),
             max(batch_table[,var], na.rm = T)))
}

get_batch_table <- function(batch_table){
    full_summary <- c()
    for(var in c('total_sum', 'NG50', 'complete', 'fragmented', 'Ns')){
        full_summary <- c(full_summary, var_summary(batch_table, var))
    }
    return(full_summary)
}
