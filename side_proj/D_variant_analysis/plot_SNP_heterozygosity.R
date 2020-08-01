library('RColorBrewer')
library('AsexStats')

tab_filename <- paste0('data/SNP_calls/', sp ,'_reduced_variants.tsv')
variant_tab <- read.table(tab_filename, stringsAsFactors = F)
colnames(variant_tab) <- c('scf', 'pos', 'qual', paste0('g', 1:5), paste0('d', 1:5))

