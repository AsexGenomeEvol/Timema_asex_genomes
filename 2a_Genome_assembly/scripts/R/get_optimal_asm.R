
# to load this file from anywhere
# source(paste0(Sys.getenv("TROOT"), '/scripts/R/get_optimal_asm.R'))

#asm <- read.table('stats/assemblies/1_Tdi_scfs.tsv', header= T)

get_asm_scores <- function(asm){
  N50_scores <- asm$N50 / max(asm$N50, na.rm=TRUE)
  NG50_scores <- asm$NG50 / max(asm$NG50, na.rm=TRUE)
  neg_size_score <- asm$diff_in_sum
  neg_size_score[neg_size_score > 0] <- 0
  neg_size_score <- 1 - abs(neg_size_score) / 430000000 # if the assembly is half of size as expected, penalise by 1
  pos_size_score <- asm$diff_in_sum
  pos_size_score[pos_size_score < 0] <- 0
  pos_size_score <- 1 - abs(pos_size_score) / 1300000000 # if the assembly is twice size as expected, penalise by 1


  asm$score <- rowMeans(cbind(N50_scores, NG50_scores, neg_size_score, pos_size_score), na.rm=TRUE)
  return(asm)
}

get_opt_asm <- function(asm){
  asm <- get_asm_scores(asm)
  return(which.max(asm$score))
}
