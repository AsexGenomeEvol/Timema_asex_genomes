library(AsexStats)

args <- commandArgs(trailingOnly=TRUE)
sp <- args[1]
if(length(args) == 1){
    batch = 'kmergenie'
} else {
    batch = args[2]
}

troot <- Sys.getenv('TROOT')
source(paste0(troot, '/scripts/R/batch_subset.R'))

if(batch == 'kmergenie' | batch == 3){
    kmergenie_file <- paste0(troot, '/stats/reads/kmertable.tsv')
    if(file.exists(kmergenie_file)){
        # if the optimal kmer is calculated then get it form computed data
        kmertab <- read.table(kmergenie_file, header = T)
        opt_k <- kmertab$opt_kmer[sp == as.character(kmertab$sp)]
    } else {
        # this data are part of the R package AsexStats
        opt_k <- timemas$optimal_kmer[timemas$codes == sp]
    }
} else if (batch == 1){
    asm <- read.table(paste0(troot, '/stats/assemblies/',sp,'_scfs.tsv'), header = T)
    opt_k <- batch_subset(asm, batch)$kmer[1]
} else if (batch == 'kmergenie.filt'){
    kmertab <- read.table(paste0(troot, '/stats/reads/kmertable_filt.tsv'), header = T)
    opt_k <- kmertab$opt_kmer[sp == as.character(kmertab$sp)]
} else {
    opt_k <- NA
}

cat(opt_k)#,"\n")
