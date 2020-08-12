library(AsexStats)
library(RColorBrewer)

# ######
# # args <- commandArgs(trailingOnly=TRUE)
# # sp_i <- 7

# asexuals <- seq(1, 10, by=2)
# to_analyse <- timemas$codes[asexuals]
sp = '5_Tge'

### LOADING SNPs
tab_filename <- paste0('data/SNP_calls/', sp, '_reduced_filtered_variants.tsv')
variant_tab <- read.table(tab_filename, stringsAsFactors = F)
colnames(variant_tab) <- c('scf', 'pos', 'qual', paste0('g', 1:5), paste0('d', 1:5))
# > head(variant_tab, 2)
#                      scf    pos    qual  g1  g2  g3  g4  g5 d1 d2 d3 d4 d5
# 1 5_Tge_b3v08_scaf000001  54809  742.15 1/1 0/0 0/0 0/0 0/0 21 19 23 20 19
# 2 5_Tge_b3v08_scaf000001 225324 3146.15 1/1 0/0 1/1 1/1 1/1 20 16 32 16 20

### LOADING INFO ABOUT GENOMES
scf_len_file <- paste0("data/", sp, "/reference/", sp, "_b3v08_scf.lengths")
scf_lengths <- read.table(scf_len_file, header = F, row.names = 1, col.names = c('scf', 'len'))
# scf_to_analyze <- lapply(scf_lengths, function(x){rownames(x)[x[,1] > 100000]})
# > head(scf_lengths)
#                            len
# 5_Tge_b3v08_scaf000001 1405450
# 5_Tge_b3v08_scaf000002 1258996
# source('../timema_assembly/scripts/prepare_reference.R')

reference <- read.table('data/external_ref/sex_lg_assigment_scores_1.4a.tsv')
colnames(reference) <- c('scf_o', 'scf', 'score', 'cov', 'len', 'asignment')
reference$chromosome <- sapply(strsplit(reference$scf, "_"), function(x) { x[1] } )

chromosomes <- data.frame(chr = paste0('lg', c(1:12, 'X')))

get_chr_size <- function(chr){
    # when scaffolds are mapped to reference, I always leave 10000 bases between them on a linkage group (consistent with the NCBI reference that used this arbitrary number of bases)
    sum(reference[reference$chromosome == chr,'len']) + ((sum(reference$chromosome == chr) - 1) * 10000)
}

chromosomes$len <- sapply( chromosomes$chr, get_chr_size )

###############################
##### plotting parameters #####
###############################

window = 1e6
gap_beween_chromosomes = 3

###############################

# I need to round up chromosome lengths to avoid one window being part of two chromosomes
chromosomes$rounded_len <- ceiling(chromosomes$len / window) * window + (window * gap_beween_chromosomes)
chromosomes$adjustments <- cumsum(c(0, chromosomes$rounded_len[1:(nrow(chromosomes) - 1)]))
rownames(chromosomes) <- chromosomes$chr
chromosomes <- chromosomes[chromosomes$chr != 'lgX', ]

get_lg_windows <- function(i) {
    all_windows <- seq(0, chromosomes[i,'rounded_len'], by = window)
    adj <- chromosomes[i, 'adjustments']
    lg_from <- all_windows[1:(length(all_windows) - 1)]
    lg_to <- all_windows[2:length(all_windows)]
    data.frame(chr = chromosomes[i, 'chr'],
               lg_from = lg_from,
               lg_to = lg_to,
               genome_from = lg_from + adj,
               genome_to = lg_from + adj)
}

variant_density_table <- do.call("rbind", lapply(1:12, get_lg_windows))
variant_density_table$variants <- 0


