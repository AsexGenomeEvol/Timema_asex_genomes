library(VennDiagram)
library(AsexStats)
source('F_structural_variation_analysis/load_SV_calls.R')
source('F_structural_variation_analysis/plot_SV_barplots.R')
source('F_structural_variation_analysis/filter_SV_calls.R')

# load SV calls
SV_call_manta_files <- paste0("data/", timemas$codes, "/variant_calls/", timemas$codes, "_survivor_manta_calls_union.vcf")
SV_calls <- load_SV_calls(SV_call_manta_files)

sv_genotype_tables <- lapply(SV_calls, get_sv_table, genotypes = T)

SV_table <- sv_genotype_tables[[1]]

het <- rowSums(SV_table[,c('00','01','02','03','04','05')] == "0/1")
SV_table$het <- het > 0

SV_table <- SV_table[SV_table$het, ]
logic_mat <- SV_table[,c('00','01','02','03','04','05')]
logic_mat[logic_mat != "./."] <- T
logic_mat[logic_mat == "./."] <- ""
logic_mat <- logic_mat[,1:5] # kill 6th individual
for( i in 1:5){
    setname <- c('A', 'B', 'C', 'D', 'F')[i]
    logic_mat[logic_mat[,i] == "TRUE", i] <- setname
}

table(apply(logic_mat, 1, function(x) { paste0(x, collapse = "") } ))

draw.quintuple.venn()