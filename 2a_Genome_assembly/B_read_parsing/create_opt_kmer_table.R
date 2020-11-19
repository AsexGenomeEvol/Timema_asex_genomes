args = commandArgs(trailingOnly=TRUE)

library(AsexStats)

kmertable_file <- 'stats/reads/kmertable'
kmergenie_output <- 'trimmed_reads/kmergenie'

if('--diploid' %in% args){
    kmertable_file <- paste0(kmertable_file, '_diploid')
    kmergenie_output <- paste0(kmergenie_output, '_diploid')
}
kmertable_file <- paste0(kmertable_file, '.tsv')
kmergenie_output <- paste0(kmergenie_output, '/histograms.dat')

# read kmergenie output for filtered reads
kmergenie <- list()
if(file.exists(kmertable_file)){
    kmer_table <- read.table(kmertable_file, header = T)
} else {
    # opt_kmer is the optimality defined by kmergenie - max genomic kmers
    kmer_table <- data.frame(sp = timemas$codes, opt_kmer = NA)
}

for(i in 1:10){
    file <- paste0('data/', timemas$codes[i], '/', kmergenie_output)
    # do not rewite already computed values
    if(is.na(kmer_table$opt_kmer[i]) & file.exists(file)){
        kmergenie[[i]] <- read.table(file, header = T)

        max_genomic_kmers <- which.max(kmergenie[[i]]$genomic.kmers)
        kmer_table$opt_kmer[i] <- kmergenie[[i]]$k[max_genomic_kmers]
    }
}

write.table(kmer_table, kmertable_file,
            quote=FALSE, sep='\t', row.names = F)
