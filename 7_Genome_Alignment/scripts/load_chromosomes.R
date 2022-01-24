get_chr_size <- function(chr){
    # when scaffolds are mapped to reference, I always leave 10000 bases between them on a linkage group (consistent with the NCBI reference that used this arbitrary number of bases)
    sum(reference[reference$chromosome == chr,'len']) + ((sum(reference$chromosome == chr) - 1) * 10000)
}

reference <- read.table('tables/sex_lg_assigment_scores_1.4a.tsv')
colnames(reference) <- c('scf_o', 'scf', 'score', 'cov', 'len', 'asignment')
reference$chromosome <- sapply(strsplit(reference$scf, "_"), function(x) { x[1] } )

chromosomes <- data.frame(chr = paste0('lg', c(1:12, 'X')))
chromosomes$len <- sapply( chromosomes$chr, get_chr_size )

# I need to round up chromosome lengths to avoid one window being part of two chromosomes
chromosomes$rounded_len <- ceiling(chromosomes$len / window) * window + (window * gap_beween_chromosomes)
chromosomes$adjustments <- cumsum(c(0, chromosomes$rounded_len[1:(nrow(chromosomes) - 1)]))
rownames(chromosomes) <- chromosomes$chr
chromosomes <- chromosomes[chromosomes$chr != 'lgX', ]
