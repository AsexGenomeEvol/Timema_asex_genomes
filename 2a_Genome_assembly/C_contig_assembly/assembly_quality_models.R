#!/usr/bin/env Rscript

# argument is species
sp = commandArgs(trailingOnly=TRUE)

# load functions, maybe I should load them individually
for(funct in dir('../scripts/R', patter = '.R', full.names=T)){
  source(funct)
}

# READ CSV of species

make_sp_model <- function(reponse, sp, asm_stats = asm_stats, model_modificatons = ''){
  model_formula <- paste0(reponse, " ~ soft + fasteris + fewdata + pse + mse + mpe + corrected", model_modificatons)
  #print(model_formula)
  model_data <- asm_stats[asm_stats$sp == sp,]
  sp_model <- lm(model_formula, data = model_data)
}


for(stat in c('NG50', 'N50')){
  assign(paste0('m_', sp, stat), make_sp_model(stat, sp, asm_stats))
}
stat <- 'diff_in_sum'
assign(paste0('m_', sp, stat), make_sp_model(stat, sp, asm_stats, ' -1'))

for(model in ls(pattern = "m_[12]")){
  print(model)
  print(summary(get(model)))
}

border <- 25
point_labels <- paste("k",asm_stats$kmer,"kc",asm_stats$kc, sep = "_")

# png(paste0('figures/',args[2],'_assemblies.png'))
#   plot(asm_stats$N50, asm_stats$NG50, xlim = c(min(asm_stats$N50) - border, max(asm_stats$N50) + border),
#        xlab = 'N50', ylab = 'NG50')
#   points(asm_stats$N50[asm_stats$corrected == 1],
#          asm_stats$NG50[asm_stats$corrected == 1],
#          pch = 20, col = 'red')
#   points(asm_stats$N50[asm_stats$fewdata == 1],
#          asm_stats$NG50[asm_stats$fewdata == 1],
#          pch = 1, col = 'blue', cex = 1.5)
#   legend('topleft', pch = c(1,20), col = c('blue','red'), legend = c('fewdata', 'corrected'))
#   text(asm_stats$N50, asm_stats$NG50, point_labels, cex = 0.8, pos = 1)
# dev.off()
