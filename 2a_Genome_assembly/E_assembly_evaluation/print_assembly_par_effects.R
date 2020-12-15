stat_type = commandArgs(trailingOnly=TRUE)[1]
sp = commandArgs(trailingOnly=TRUE)[2]
#print(stat_type)
#print(sp)

library(AsexStats)

# takes optimal kmer assembly only
collapse_kmer_kc_par <- function(sp_stats, parameters, aggregration_var = 'NG50', funct = max, ...){
  out <- tryCatch(
    {
      fomula <- paste0(aggregration_var, ' ~ ', paste0(parameters, collapse = ' + '))
      aggregration_var_values <- aggregate(as.formula(fomula), FUN=funct, data=sp_stats, ...)[,aggregration_var]
      sp_stats_reduced <- sp_stats[sp_stats[,aggregration_var] %in% aggregration_var_values,]
      return(sp_stats_reduced)
    },
    error=function(cond) {
      return(sp_stats)
    }
  )
  return(out)
}

# find the sole effects of parameter change
get_par_improvements <- function(sp_stats){
  # col sp : SP, NA, NA or SP, SP, SP
  # col rows: N50, NG50, Total
  # cols:
  parameters <- c("soft", "fasteris", "fewdata", "mse", "pse", "mpe", "corrected")
  variables <- c('sp', 'metric', parameters)
  metrics <- c('N50', 'NG50', 'diff_in_sum')

  absent_vars <- c(parameters, metrics)[!(c(parameters, metrics) %in% colnames(sp_stats))]
  if(length(absent_vars) != 0){
    print('I miss folowing variables:')
    print(absent_vars)
    stop()
  }

  sp_stats <- collapse_kmer_kc_par(sp_stats, parameters, 'NG50', max, rm.na = T)
  sp_improvements <- make_data_frame(variables)

  sp_improvements[1:3,'sp'] <- rep(sp_stats$sp[1], 3)
  sp_improvements[1:3,'metric'] <- metrics
  for(parameter in parameters){
    # select parameters for grouping (all but not the one that is evaluated)
    other_parameters <- parameters[parameters != parameter]
    # make a stirng of the status of the grouping parameters (easier for comparison)
    sp_stats$grouping <- apply(sp_stats[,other_parameters], 1,
                               FUN = function(x){paste0(x, collapse = '')})
    group_table <- table(sp_stats$grouping)
    group_table <- group_table[group_table >= 2]
    stat_diff <- data.frame()

    for(group in names(group_table)){
      group_stats <- sp_stats[sp_stats$grouping == group, ]
      alt_values <- as.numeric(group_stats[,parameter])
      # for any T/F value
      if(!any(alt_values == 1)){
        next
      }
      alt <- 0
      # for soft
      if(any(alt_values >= 2)) alt <- max(alt_values)
      stat_diff <- rbind(stat_diff,
                         abs(group_stats[as.numeric(group_stats[,parameter]) == 1, metrics]) -
                         abs(group_stats[as.numeric(group_stats[,parameter]) == alt, metrics]))
    }

    if(nrow(stat_diff) > 0){
      chosen_comparison <- which.max(stat_diff$NG50)
      sp_improvements[1:3,parameter] <- as.numeric(stat_diff[chosen_comparison,])
    }
  }
  return(sp_improvements)
}

files_of_type <- dir('assemblies', pattern = stat_type, full.names=T)
sp_file <- files_of_type[grep(sp, files_of_type)]
sp_stats <- read.table(sp_file, header = T)
sp_stats$sp <- sp
print(sp)
print(get_par_improvements(sp_stats))
