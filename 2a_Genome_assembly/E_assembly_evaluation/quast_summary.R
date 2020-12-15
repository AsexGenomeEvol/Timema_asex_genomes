# wild script (not called by makefile)

source('scripts/R/variables.R')

quast_files <- c(
"data/1_Tdi/assembly/abyss_k87_nomse_nompe_fewdata_BESST/Scaffolds_pass5_GC_quast/report.tsv",
"data/1_Tps/assembly/abyss_k83_nomse_nompe_fewdata_BESST/Scaffolds_pass5_GC_quast/report.tsv",
"data/2_Tsi/assembly/abyss_k87_nomse_nompe_fewdata_BESST/Scaffolds_pass5_GC_quast/report.tsv",
"data/2_Tcm/assembly/abyss_k83_nomse_nompe_fewdata_BESST/Scaffolds_pass5_GC_quast/report.tsv",
"data/3_Tms/assembly/abyss_k89_nomse_nompe_fewdata_BESST/Scaffolds_pass5_GC_quast/report.tsv",
"data/3_Tce/assembly/abyss_k83_nomse_nompe_fewdata_BESST/Scaffolds_pass5_GC_quast/report.tsv",
"data/4_Tte/assembly/abyss_k81_nomse_nompe_fewdata_BESST/Scaffolds_pass5_GC_quast/report.tsv",
"data/4_Tbi/assembly/abyss_k81_nomse_nompe_fewdata_BESST/Scaffolds_pass5_GC_quast/report.tsv",
"data/5_Tge/assembly/abyss_k87_nomse_nompe_fewdata_BESST/Scaffolds_pass5_GC_quast/report.tsv",
"data/5_Tpa/assembly/abyss_k65_nomse_nompe_fewdata_BESST/Scaffolds_pass5_GC_quast/report.tsv")

quast_list <- list()
for(i in 1:10){
    quast_list[[i]] <- read.table(header = T, quast_files[i], comment.char='!', sep = '\t', quote="")
}

plot_barplot <- function(vector, title){
    barplot(vector, main = title, col = c(asex_blue, sex_red), width = 0.9)
    text(seq(0.25,10.5, length = 10), par("usr")[3] - 5,
         srt = 20, pos = 1, xpd = TRUE, col = c(asex_blue, sex_red),
         labels = timema_labels)
    legend('topright', col = c(asex_blue, sex_red), pch = 20, c('asex','sex'), bty = 'n')
}

pdf('stats/figures/b3v04/asm_length.pdf')
    plot_barplot(round((unlist(lapply(quast_list, function(x){x[10,2]}))) / 1e6, 1), 'sum of scaffolds (>10kbp)')
    mtext('sum of scfs [Mbp]', side = 2, line = +2)
dev.off()

pdf('stats/figures/b3v04/predicted_long_genes.pdf')
    plot_barplot(unlist(lapply(quast_list, function(x){x[26,2]})), 'long predicted genes (>3kbp)')
dev.off()

pdf('stats/figures/b3v04/predicted_all_genes.pdf')
    plot_barplot(unlist(lapply(quast_list, function(x){x[24,2]})), 'all predicted genes (>300bp)')
dev.off()
