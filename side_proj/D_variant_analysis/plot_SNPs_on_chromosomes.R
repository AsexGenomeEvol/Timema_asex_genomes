library(AsexStats)
library(RColorBrewer)

asexuals <- seq(1, 10, by=2)
to_analyse <- timemas$codes[asexuals]

# load genome mapping to cristinae reference
aln_files <- paste0('data/b3v08_anchoring_to_LGs/', to_analyse, '_scf_block_alignment.tsv')
genome_aligments <- lapply(aln_files, read.table, header = T, stringsAsFactors = F)

### LOADING INFO ABOUT GENOMES
scf_len_files <- paste0("data/", to_analyse, "/reference/", to_analyse, "_b3v08_scf.lengths")
scf_lengths <- lapply(scf_len_files, read.table, header = F, row.names = 1, col.names = c('scf', 'len'))
scf_to_analyze <- lapply(scf_lengths, function(x){rownames(x)[x[,1] > 100000]})

### LOADING SNPs
hetero_SNP_files <- paste0('Timema_SNP_calling/data/', to_analyse, '_heterozygous_SNP_filter_passed.tsv')
homo_SNP_files <- paste0('Timema_SNP_calling/data/', to_analyse, '_homozygous_SNP_filter_passed.tsv')

heteroz_variants <- lapply(hetero_SNP_files, read.table, col.names=c('scf', 'pos'), stringsAsFactors = F)
homo_variants <- lapply(homo_SNP_files, read.table, col.names=c('scf', 'pos'), stringsAsFactors = F)

# ######
# # args <- commandArgs(trailingOnly=TRUE)
# # sp_i <- 7
source('../timema_assembly/scripts/prepare_reference.R')

