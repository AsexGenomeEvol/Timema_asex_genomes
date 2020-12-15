library(AsexStats)
stats <- read.table('stats/reference/genome_table.tsv', header = T, sep = "\t")



#### Asm length

pdf('figures/b3v07/asm_sizes.pdf')

par(mar=c(c(5, 5, 1.2, 1) + 0.1))
barplot(stats$total_sum, col = c(asex_blue, sex_red), width = 0.9, ylim = c(0, 1.381 + 0.05),
        names.arg = NA, cex.axis = 1.3, cex.lab = 2, ylab = "Assembly size [Gbp]")
text(seq(0.3,10.2, length = 10) - 0.5, par("usr")[3] - 0.08,
     srt = 30, pos = 1, xpd = TRUE, col = c(asex_blue, sex_red),
     labels = timemas$labels, cex = 1.4)
lines(c(0,11), c(1.381, 1.381), lwd = 2)
text(8.6, 1.36, labels = "flow cytometry estimate", pos = 3, cex = 1.4)

dev.off()

#### N50

pdf('figures/b3v07/N50.pdf')

par(mar=c(c(5, 5, 1.2, 1) + 0.1))
barplot(stats$N50, col = c(asex_blue, sex_red), width = 0.9,
        names.arg = NA, cex.axis = 1.3, cex.lab = 2, ylab = "N50 [kbp]")
text(seq(0.3,10.2, length = 10) - 0.5, par("usr")[3] - 8,
     srt = 30, pos = 1, xpd = TRUE, col = c(asex_blue, sex_red),
     labels = timemas$labels, cex = 1.3)


dev.off()

#### BUSCO (stacked barplot?)

pdf('figures/b3v07/BUSCO.pdf')

light_asex_blue <- "#92C5DE8D"
light_sex_red <- "#D6604DA2"
dark_asex_blue <- "#92C5DEFF"
dark_sex_red <- "#D6604DFF"

par(mar=c(c(5, 5, 1.2, 1) + 0.1))
barplot(stats$complete + stats$fragmented, col = c(light_asex_blue, light_sex_red), width = 0.9,
        names.arg = NA, cex.axis = 1.2, cex.lab = 2, ylab = "BUSCO [%]", ylim = c(0,100))
barplot(stats$complete, col = c(asex_blue, sex_red), width = 0.9, add = T, yaxt="n")
barplot(stats$duplicated, col = c(dark_asex_blue, dark_sex_red), width = 0.9, add = T, yaxt="n")
text(seq(0.3,10.2, length = 10) - 0.5, par("usr")[3] - 6,
     srt = 30, pos = 1, xpd = TRUE, col = c(asex_blue, sex_red),
     labels = timemas$labels, cex = 1.3)
# text(par("usr")[2] - 1, c(45, 72),
#      srt = 90, pos = 1, xpd = TRUE,
#      labels = c("complete", "fragmented"), cex = 1.8)
text(0.5, c(88, 94.5), cex = 1.8, labels = c('C','F'))

dev.off()
