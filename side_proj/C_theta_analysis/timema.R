##### READ_DATA
args <- commandArgs(trailingOnly=TRUE)

library(AsexStats)

ind <- args[1]
ref <- args[2]
window <- args[3]
ind <- 'ref'
ref <- 'b3v05'
window <- '100000'

logThetas <- list()
atlasData <- list()

# PREPARED FOR MORE ESTIMATES
# for(i in seq(3, 4, by = 2)){
#     asex_sp = timemas[i]
#     sex_sp = timemas[i+1]
#
#     pdf(paste0('variant_stats/figures/filt_loghist_', asex_sp, '_', sex_sp, '.pdf'))
#
#     sex_data <- load_thetas(i)
#     asex_data <- load_thetas(i+1)
#
#     logThetas[[i]] <- log10(asex_data$theta_MLE)
#     logThetas[[i+1]] <- log10(sex_data$theta_MLE)
#
#     plot_log_hist(sex_data = sex_data$theta_MLE, asex_data = asex_data$theta_MLE,
#                   breaks = 160, barwidth = 4, xlab = expression(theta))
#     legend('topright', c(sex_sp, asex_sp),
#            col = c(sex_red, asex_blue), pch = 20, bty = 'n')
#
#     dev.off()
# }

for(i in 1:10){
    # non-filtered all data
    atlasData[[i]] <- load_thetas(timemas$codes[i], ind, ref, window)
    # filtered theta estimates
    logThetas[[i]] <- log10(filter_theta(atlasData[[i]], filt_cov = F, window_size = 9999)[,'theta_MLE'])
}

#require(vioplot)
#require(digest)
library(sm)
source('/Volumes/dump/projects/timema/AsexStats/almost_working/vioplot2.R')
#source('scripts/external/vioplot2.R')
theta_label <- expression(paste("Heterozygosity of 100kbp windows [" , log[10], " " , theta, ']'))

quartz(type = 'pdf', file = 'figures/tim_w100000_theta_violin.pdf')
plot(x=NULL, y=NULL,
     xlim = c(0.5, 5.5), ylim=c(min(unlist(logThetas)), max(unlist(logThetas))),
     type="n", ann=FALSE, axes=F, bty="n")
axis(1, at=c(1:5), labels = F)
text(c(1:5), par("usr")[3] - 0.2, srt = 10, pos = 1, xpd = TRUE,
     labels = species_pairs$labels)
axis(2)
mtext(theta_label, side = 2, line = +2, cex = 1.3)

for(i in c(7:8,5:6,1:4,9:10)){
    vioplot2(logThetas[[i]],
             at = ifelse(i %% 2 == 1, (i + 1) / 2, i / 2),
             side = ifelse(i %% 2 == 1, "left", "right"),
             col = ifelse(i %% 2 == 1, asex_blue, sex_red),
             add = T, h = 0.2)
}
# sex_legend(cex = 1.3)
dev.off()

quartz(type = 'pdf', file = 'figures/tim_w100000_theta_violin_zoomed.pdf')
plot(x=NULL, y=NULL,
     xlim = c(0.5, 5.5), ylim=c(-4, max(unlist(logThetas))),
     type="n", ann=FALSE, axes=F, bty="n")
axis(1, at=c(1:5), labels = F)
text(c(1:5), par("usr")[3] - 0.2, srt = 10, pos = 1, xpd = TRUE,
     labels = species_pairs$labels)
axis(2)
mtext(theta_label, side = 2, line = +2, cex = 1.3)

for(i in c(7:8,5:6,1:4,9:10)){
    to_plot <- logThetas[[i]][logThetas[[i]] > -4]
    vioplot2(to_plot,
             at = ifelse(i %% 2 == 1, (i + 1) / 2, i / 2),
             side = ifelse(i %% 2 == 1, "left", "right"),
             col = ifelse(i %% 2 == 1, asex_blue, sex_red),
             add = T, h = 0.2)
}

# sex_legend(cex = 1.3)
dev.off()
asex_indices <- c(7:8,5:6,1:4,9:10)[seq(1,10,by=2)]

medians <- list()
means <- list()

for(i in 1:10){
    medians[[i]] <- median(10^(logThetas[[i]]))
    means[[i]] <- mean(10^(logThetas[[i]]))
}

asex_medians_pair_ord <- unlist(medians)[asex_indices]
asex_means_pair_ord <- unlist(medians)[asex_indices]

#species_pairs$divergence_times

#plot(unlist(means[seq(1,10,by=2)]) ~ divergence_times, pch = 20)
# unlist(lapply(logThetas[seq(1,10,by=2)], function(x){ sum(x < -5)})) / 1000
pdf('figures/heterozygosity_mean_vs_age.pdf')
par(mar=c(5.1,4.2,0,2))
plot(asex_means_pair_ord ~ species_pairs$divergence_times,
     pch = 20, xlab = 'Mean species pair divergence [%]', ylab = 'mean theta [2Tmu]',
     cex.lab = 1.3, cex = 1.5, bty = 'n', xlim = c(0.7, 1.95)) #, , ylim = c(15, 21.1)
     text(asex_means_pair_ord ~ species_pairs$divergence_times,
         labels = timemas$labels[asex_indices], pos = 3)
dev.off()

pdf('figures/heterozygosity_median_vs_age.pdf')
par(mar=c(5.1,4.2,0,2))
plot(asex_medians_pair_ord ~ species_pairs$divergence_times,
     pch = 20, xlab = 'Mean species pair divergence [%]', ylab = 'median theta [2Tmu]',
     cex.lab = 1.3, cex = 1.5, bty = 'n', xlim = c(0.7, 1.95)) #, ylim = c(15, 21.1)
     text(asex_medians_pair_ord ~ species_pairs$divergence_times,
          labels = timemas$labels[asex_indices], pos = 3)
dev.off()



# #2.8 × 10(-9) / site / generation
# 300*10^6
# 2.8 * 1.3 * 2 * 10e6
# # 7 280 000
# # 0.00028 after 0.5 * 10e6
# # 0.00056 after 1 * 10e6
# # 0.00112 after 2 * 10e6
# 0.00095 - 0.00120
