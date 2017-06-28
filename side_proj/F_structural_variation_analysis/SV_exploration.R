#library('vcfR')

args <- commandArgs(trailingOnly=TRUE)
ind <- args[1]
ref <- args[2]
# ind <- 'ref_is350'
# ref <- 'b3v04'

source('scripts/R/variables.R')
source('scripts/R/ssplit.R')
source('scripts/R/sex_legend.R')
# source('scripts/R/plot_log_hist.R')
source('variant_stats/load_sv.R')

ind <- 'ref_is350'
ref <- 'b3v04'

SV_colnames <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE1')

# sex_sp <- '4_Tbi'
# sex_file <- paste0('data/', sex_sp , '/variant_calls/ref_is350/manta/results/variants/diploidSV.vcf')
# sex_vcf <- read.vcfR( sex_file, verbose = FALSE )

# note, capital letters just respect original naming conventions of the VCF file

sv_list <- list()
for(i in c(1:10)){
    sv_list[[i]] <- load_sv(timemas[i])
}

sp_tables <- list()
for(i in c(1:10)){
    sp_tables[[i]] <- table(sv_list[[i]]$SVTYPE)[-1]
}

all_tables <- unlist(sp_tables)

# Grouped Bar Plot
pdf('variant_stats/figures/ref_is350_barplot.pdf')
par(mfrow = c(2,5))
ylim <- c(0, max(unlist(all_tables)))
for(i in c(1,3,5,7,9)){
    barplot(sp_tables[[i]], col = asex_blue, main = timema_labels[i], ylim = ylim, cex.axis = 1.6, cex.names = 1.6)
}
for(i in c(2,4,6,8,10)){
    barplot(sp_tables[[i]], col = sex_red, main = timema_labels[i], ylim = ylim, cex.axis = 1.6, cex.names = 1.6)
}
dev.off()

for(i in 1:4){
    variant <- unique(names(all_tables))[i]
    pdf(paste0('variant_stats/figures/', ind, '_', ref, '_', variant, '_barplot.pdf'))
#    variant <- 'DEL'
        variant_label <- c('Deletions', 'Duplications', 'Insertions', 'Inversions')[i]
    one_table <- all_tables[names(all_tables) == variant]
    names(one_table) <- timemas[1:10]
    barplot(one_table, main = variant_label, col = c(asex_blue, sex_red), width = 0.9,
            names.arg = NA, cex.axis = 1.2, cex.main = 2)
    text(seq(0.2,10.45, length = 10), par("usr")[3] - (max(one_table) * 0.02),
         srt = 20, pos = 1, xpd = TRUE, col = c(asex_blue, sex_red),
         labels = timema_labels, cex = 1.3)
#    legend('topright', col = c(asex_blue, sex_red), pch = 20, c('asex','sex'), bty = 'n')
    dev.off()
}


## Relate to SV length
del_sizes <- list()
dup_sizes <- list()
ins_sizes <- list()
inv_sizes <- list()
#DEL  DUP  INS  INV

for(i in c(1:10)){  # already computed timemas
    del_sizes[[i]] <- sv_list[[i]][sv_list[[i]]$SVTYPE == 'DEL','SVLEN']
    dup_sizes[[i]] <- sv_list[[i]][sv_list[[i]]$SVTYPE == 'DUP','SVLEN']
    ins_sizes[[i]] <- sv_list[[i]][sv_list[[i]]$SVTYPE == 'INS','SVLEN']
    inv_sizes[[i]] <- sv_list[[i]][sv_list[[i]]$SVTYPE == 'INV','SVLEN']
}

# sizes <- del_sizes
# svtype <- 'DEL'

require(vioplot)
require(digest)
library(sm)
source('scripts/external/vioplot2.R')

############
#
###########
plotLogViolins <- function(sizes){
    ylim <- c(min(log10(abs(unlist(sizes))), na.rm = T),
              max(log10(abs(unlist(sizes))), na.rm = T))
    plot(x=NULL, y=NULL,
         xlim = c(0.5, 5.5), ylim = ylim,
         type="n", ann=FALSE, axes=F)
    axis(1, at=c(1:5), labels = F)
    text(c(1:5), par("usr")[3] - 0.2, srt = 10, pos = 1, xpd = TRUE,
         labels = c(expression(italic("T. douglasi / T. poppensis")),
                    expression(italic("T. shepardi / T. californicum")),
                    expression(italic("T. monikensis / T. cristinae")),
                    expression(italic("T. tahoe / T. barmani")),
                    expression(italic("T. genevieve / T. podura"))))
    mtext(expression(paste("Structural variant sizes [" , log[10], " nt]")), side = 2, line = +2, cex = 1.3)


    for(i in 1:10){
        vioplot2(log10(abs(sizes[[i]])),
                 at = ifelse(i %% 2 == 1, (i + 1) / 2, i / 2),
                 side = ifelse(i %% 2 == 1, "left", "right"),
                 col = ifelse(i %% 2 == 1, asex_blue, sex_red),
                 add = T, h = 0.2)
    }
}

pdf(paste0('variant_stats/figures/', ind, '_', ref, '_DEL_violin_sizes.pdf'))
    plotLogViolins(del_sizes)
    axis(2, at = c(1, 2, 3, 4, 5), labels = c('10', '100', '1000', '10000', '100000'), cex.axis = 1)
    title('Deletions')
dev.off()

pvals <- c()
#tpvals <- c()
for(i in seq(1,10,by=2)){
    pvals <- c(pvals, wilcox.test(dup_sizes[[i]], dup_sizes[[i+1]])$p.value)
#    tpvals <- c(tpvals, t.test(dup_sizes[[i]], dup_sizes[[i+1]])$p.value)
}

pdf(paste0('variant_stats/figures/', ind, '_', ref, '_DUP_violin_sizes.pdf'))
    plotLogViolins(dup_sizes)
    axis(2, at = c(1, 2, 3, 4, 5), labels = c('10', '100', '1000', '10000', '100000'), cex.axis = 1)
    title('Duplications')
    # text(0.5, 5.4, labels = 'pval')
    # text(1:5, 5.4, labels = format(pvals, digits=1))
dev.off()


plotLogSizes <- function(sizes, svtype){
    variant_label <- c('Deletions', 'Duplications', 'Insertions', 'Inversions')[svtype == c('DEL','DUP','INS','INV')]
    ylim <- c(min(log10(abs(unlist(sizes))), na.rm = T),
              max(log10(abs(unlist(sizes))), na.rm = T))
    plot(numeric(0),
         xlim = c(0.5,10.5), ylim = ylim, ann=FALSE, axes=F)
    axis(1, at=c(1:10), labels = F)
    text(c(1:10), par("usr")[3] - ((ylim[2] - ylim[1]) / 20), srt = 15, pos = 1, xpd = TRUE, col = c(asex_blue, sex_red),
         labels = timema_labels)
    axis(2)
    mtext(expression(paste(log[10], ' structural variant sizes')), side = 2, line = +2, cex = 1.3)
    title(variant_label)
    for(i in c(1:10)){
        col <- ifelse(i %% 2 == 1, asex_blue, sex_red)
        boxplot(log10(abs(sizes[[i]])), at = i, add = T, pch = 20, col = col)
    }
}
#sex_legend
pdf(paste0('variant_stats/figures/', ind, '_', ref, '_DEL_sizes.pdf'))
    plotLogSizes(del_sizes, 'DEL')
dev.off()

pdf(paste0('variant_stats/figures/', ind, '_', ref, '_INV_sizes.pdf'))
    plotLogSizes(inv_sizes, 'INV')
dev.off()

pdf(paste0('variant_stats/figures/', ind, '_', ref, '_DUP_sizes.pdf'))
    plotLogSizes(dup_sizes, 'DUP')
dev.off()

pdf(paste0('variant_stats/figures/', ind, '_', ref, '_INS_sizes.pdf'))
    plotLogSizes(ins_sizes, 'INS')
dev.off()


# axis timemas[1:8]

## show no relation to quality of assembly
# no correlation between quality of genome and # SVs
genome_summary <- read.table('stats/reference/b3v04.tsv', header = T)
#plot(c(unlist(lapply(sp_tables, sum))) ~ genome_summary$complete[1:8])
plotSexAsex <- function(x, y, main = '', xlab = '', ylab = ''){
    asex <- seq(1,10, by = 2)
    sex <- seq(2,10, by = 2)
    plot(x[asex], y[asex],
        xlim = c(min(x), max(x)),
        ylim = c(min(y), max(y)),
        col = asex_blue,
        pch = 20,
        cex = 1.5,
        main = main,
        xlab = xlab, ylab = ylab,
        cex.main = 2, cex.lab = 1.3
    )
    points(x[sex], y[sex],
        col = sex_red,
        pch = 20,
        cex = 1.3)
}


pdf(paste0('variant_stats/figures/', ind, '_', ref, '_SVs_vs_NG50.pdf'))
    plotSexAsex(genome_summary$NG50, c(unlist(lapply(sp_tables, sum))),
        xlab = 'NG50', ylab = 'Total structural variants')
dev.off()

pdf(paste0('variant_stats/figures/', ind, '_', ref, '_duplications_vs_NG50.pdf'))
    plotSexAsex(genome_summary$NG50, unlist(lapply(dup_sizes, length)),
        xlab = 'NG50', ylab = 'Duplications')
dev.off()


source('/Volumes/dump/scripts/R/table_print.R')

sp_sv_tab <- cbind(data.frame(sp = paste0('\\textit{', timema_names, '} ', c('\\Female', '\\Female\\Male')),
                              as.data.frame(t(matrix(nrow = 4, unlist(sp_tables))))))
colnames(sp_sv_tab) <- c('species', 'deletions', 'duplications', 'insertions', 'inversions')


print.latextable(sp_sv_tab)
