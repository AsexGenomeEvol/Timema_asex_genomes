# wild script

ref <- commandArgs(trailingOnly=TRUE)
# ref <- 'b3v04'

source('scripts/R/variables.R')
source('scripts/R/plot_log_hist.R')

for (i in seq(1,10, by = 2)) {

    asex_sp <- timemas[i]
    sex_sp <- timemas[i + 1]
    #res = 300
    pdf(paste0('stats/reference/', asex_sp, '_', sex_sp, '_', ref ,'_scf.pdf'))

    asex_file <- paste0('data/', asex_sp, '/reference/',
                        asex_sp, '_', ref, '_lengths.tsv')
    sex_file <- paste0('data/', sex_sp, '/reference/',
                        sex_sp, '_', ref, '_lengths.tsv')

    asex_lengths <- read.table(asex_file, header = F)$V1
    sex_lengths <- read.table(sex_file, header = F)$V1

    plot_log_hist(asex_lengths, sex_lengths, noborder = T,
        col = c("#92C5DECD", "#D6604DED"), barwidth = 8,
        axes = F, xlim = c(0, 2000000))
    axis(1, at = c(100000, 500000, 1000000, 2000000), labels = c('100000', '500000', '1000000', '2000000'), cex.axis = 1.6)
    axis(2, cex.axis = 1.6)

    dev.off()
}


# read.table('stats/reference/b3v04.tsv', header = T)
# print latex table
