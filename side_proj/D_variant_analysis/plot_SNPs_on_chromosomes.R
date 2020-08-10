# library(AsexStats)
library(RColorBrewer)

# asexuals <- seq(1, 10, by=2)
# to_analyse <- timemas$codes[asexuals]
sp = '5_Tge'

tab_filename <- paste0('data/SNP_calls/', sp, '_reduced_filtered_variants.tsv')
variant_tab <- read.table(tab_filename, stringsAsFactors = F)
colnames(variant_tab) <- c('scf', 'pos', 'qual', paste0('g', 1:5), paste0('d', 1:5))

# load genome mapping to cristinae reference
aln_file <- paste0('data/b3v08_anchoring_to_LGs/', sp, '_scf_block_alignment.tsv')
genome_aligment <- read.table(aln_file, header = T, stringsAsFactors = F)

### LOADING INFO ABOUT GENOMES
scf_len_file <- paste0("data/", sp, "/reference/", sp, "_b3v08_scf.lengths")
scf_lengths <- read.table(scf_len_file, header = F, row.names = 1, col.names = c('scf', 'len'))
# scf_to_analyze <- lapply(scf_lengths, function(x){rownames(x)[x[,1] > 100000]})

### LOADING SNPs
# data/<sp>/reference/<sp>_b3v08_scf.lengths
# data/b3v08_anchoring_to_LGs/<sp>_scf_block_alignment.tsv
# data/SNP_calls/<sp>_reduced_filtered_variants.tsv

# ######
# # args <- commandArgs(trailingOnly=TRUE)
# # sp_i <- 7
source('../timema_assembly/scripts/prepare_reference.R')

