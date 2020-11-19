outfile = commandArgs(trailingOnly=TRUE)[1]

library(AsexStats)

#data/1_Tdi/trimmed_reads/mp_fasteris_FR/1_Tdi_R1t_is_5000.fq.gz	480818350
# 1_Tdi_trimmed_reads.tsv

read_list_files <- paste0('stats/reads/', timemas$codes, '_trimmed_reads.tsv')

pe_libs <- paste0('is_', c('350','350_run2','550','700'))
mp_libs <- paste0('is_', c('3000_fasteris','5000_fasteris',
                          '3000_nxtrim',  '5000_nxtrim'))

coverages <- make_data_frame(c('sp', pe_libs, mp_libs,
                               'pse', 'pse_run2', 'mpe','mse'))

pe_srearch_strings <- paste0('t_', pe_libs, '.fq')

get_cov <- function(pattern, sp_table){
  sum(sp_table[grepl(pattern, sp_table$V1),2]) / 1300000000
}

for(i in 1:10){
  sp <- timemas$codes[i]
  file_to_parse <- read_list_files[i]
  sp_table <- read.table(file_to_parse)
  pse <- get_cov('np_is_....fq', sp_table)
  pse_run2 <- get_cov('np_is_350_run2', sp_table)
  mpe <- get_cov('225', sp_table)
  mse <- get_cov('_se_mp', sp_table)
  pe_cov <- c()
  for(pe in 1:4){
    pe_cov[pe] <- get_cov(pe_srearch_strings[pe], sp_table)
  }
  me_cov <- c()
  for(mptype in c('fasteris', 'nxtrim')){
    sp_subtable <- sp_table[grepl(mptype, sp_table$V1),]
    me_cov <- c(me_cov, get_cov('3000', sp_subtable))
    me_cov <- c(me_cov, get_cov('5000', sp_subtable))
  }
  coverages[i,'sp'] <- sp
  coverages[i,2:13] <- round(c(pe_cov, me_cov, pse, pse_run2, mpe, mse), 2)
}

coverages$total_pe_fewdata <- rowSums(coverages[,pe_libs[c(1,3,4)]])
coverages$total_pe <- rowSums(coverages[,pe_libs])
coverages$total_pse <- rowSums(coverages[,c('pse', 'pse_run2')])
coverages$total <- rowSums(coverages[,2:13])

write.table(coverages, outfile, quote = F, sep = '\t', row.names = F)
