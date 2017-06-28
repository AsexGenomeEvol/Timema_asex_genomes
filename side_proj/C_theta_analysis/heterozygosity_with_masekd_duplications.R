##### READ_DATA
args <- commandArgs(trailingOnly=TRUE)

ind <- args[1]
ref <- args[2]
# ind <- 'ref_is350'
# ref <- 'b3v04'

source('scripts/R/variables.R')
source('scripts/R/plot_log_hist.R')
source('scripts/R/sex_legend.R')
source('variant_stats/load_thetas.R')
source('variant_stats/load_sv.R')

# # # # # # # #
# LOAD THETA  #
# # # # # # # #

computed <- c(1,3,4,5,7,9)

atlasData <- list()
for(i in computed){
    atlasData[[i]] <- filter_theta(load_thetas(i))
}

# # # # # # #
# LOAD SVs  #
# # # # # # #

sv_list <- list()
for(i in c(1:10)){
    sv_list[[i]] <- load_sv(timemas[i])
    sv_list[[i]] <- sv_list[[i]][sv_list[[i]]$SVTYPE != 'BND',]
}

# # # # # # # # # # # #
# FILTER THETA BY SVs #
# # # # # # # # # # # #

removed_data <- list()
filtered_data <- list()
i <- 1

sv_index <- 1
for(scf in levels(atlasData[[i]]$Chr)){
    removed <- data.frame()
    kept <- data.frame()
    ch_data <- atlasData[[i]][atlasData[[i]]$Chr == scf,]
    ch_sv_data <- sv_list[[i]][sv_list[[i]]$CHROM == scf,]

    for(window in ch_data){
        if(any(window$end > ch_sv_data$POS & window$end < ch_sv_data$POS + ch_sv_data$SVLEN)){
            removed <- rbind(removed, window)
        } else {
            kept <- rbind(kept, window)
        }
    }

    removed_data[[i]] <- removed
    filtered_data[[i]] <- kept
}
