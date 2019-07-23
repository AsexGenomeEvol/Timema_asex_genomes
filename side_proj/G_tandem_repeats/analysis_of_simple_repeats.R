
library(AsexStats)

dirs <- paste('data', timemas$code, 'tandem_repeats', sep = '/')

rep_summary <- data.frame(files = unlist(lapply(dirs, dir, pattern = '_kseek.rep.total', full.names = T)), stringsAsFactors = F)
rep_annot <- lapply(rep_summary$files, read.table, col.names = c('pattern', 'rep'), stringsAsFactors = F)

rep_summary$sample <- substr(rep_summary$files, 27, 32)
rep_summary$annotated_rep <- sapply(rep_annot, nrow)
rep_summary$annotated_nt <- sapply(rep_annot, function(x) { sum(nchar(x$pattern) * x$rep) } ) / 1e6

all_motifs <- unlist(lapply(rep_annot, function(x) { x$pattern }))

motif_tab <- data.frame(motif = unique(all_motifs))
row.names(motif_tab) <- motif_tab$motif

motif_tab[, rep_summary$sample] <- 0
for ( i in 1:nrow(rep_summary) ){
    motif_tab[rep_annot[[i]]$pattern, rep_summary$sample[i]] <- rep_annot[[i]]$rep
}

# I am here excluding mononucleotide repeats as they are probably just a noise
mono_nt_tab <- subset(motif_tab, motif %in% c('A','G','C','T'))
motif_tab <- subset(motif_tab, ! motif %in% c('A','G','C','T'))

write.table(motif_tab, 'data/tandem_repeats.tsv', quote = F, sep = '\t', row.names = F)


motif_tab$found_in <- rowSums(motif_tab[,2:61] > 0)

conserved_motif_tab <- motif_tab[motif_tab$found_in == 60, ]

# nrow(conserved_motif_tab)
# 439

individual_distance <- matrix(rep(NA, 60 * 60), nrow = 60)
for ( i in 1:60 ){
    for ( j in 1:60){
        individual_distance[i, j] <- sum(motif_tab[,i+1] > 0 & motif_tab[,j+1] > 0)
    }
}

library(gplots)

colnames(individual_distance) <- rep_summary$sample
rownames(individual_distance) <- rep_summary$sample
heatmap.2(log10(individual_distance), trace = 'none')

conserved_rep_mat <- as.matrix(conserved_motif_tab[,2:61])
colnames(conserved_rep_mat) <- rep_summary$sample

heatmap.2(conserved_rep_mat, trace = 'none')
heatmap.2(log10(conserved_rep_mat), trace = 'none')