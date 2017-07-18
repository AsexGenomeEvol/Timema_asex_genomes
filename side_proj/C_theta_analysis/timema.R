##### READ_DATA
args <- commandArgs(trailingOnly=TRUE)

library(AsexStats)

# ind <- args[1]
# ref <- args[2]
# window <- args[3]
ind <- 'ref'
ref <- 'b3v05'
window <- '100000'

files <- paste0('data/',timemas$codes,'/variant_calls/ref/atlas/',
                ind, '_to_', ref, '_w', window, '_theta_estimates.txt')

atlas_data <- list()
for(i in 1:10){
    # non-filtered all data
    atlas_data[[i]] <- read.table(files[i], header = T)
}

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
#     log_thetas[[i]] <- log10(asex_data$theta_MLE)
#     log_thetas[[i+1]] <- log10(sex_data$theta_MLE)
#
#     plot_log_hist(sex_data = sex_data$theta_MLE, asex_data = asex_data$theta_MLE,
#                   breaks = 160, barwidth = 4, xlab = expression(theta))
#     legend('topright', c(sex_sp, asex_sp),
#            col = c(sex_red, asex_blue), pch = 20, bty = 'n')
#
#     dev.off()
# }

thetas <- lapply(atlas_data, function(x){filter_thetas(x, filt_cov = F, window_size = 9999)[,'theta_MLE']} )
log_thetas <- lapply(thetas, log10)

theta_label <- expression(paste("Heterozygosity of 100kbp windows [" , log[10], " " , theta, ']'))

quartz(type = 'pdf', file = 'figures/tim_w100000_theta_violin.pdf')
plot(x=NULL, y=NULL,
     xlim = c(0.5, 5.5), ylim=c(min(unlist(log_thetas)), max(unlist(log_thetas))),
     type="n", ann=FALSE, axes=F, bty="n")
axis(1, at=c(1:5), labels = F)
text(c(1:5), par("usr")[3] - 0.2, srt = 10, pos = 1, xpd = TRUE,
     labels = timema_pairs$labels)
axis(2)
mtext(theta_label, side = 2, line = +2, cex = 1.3)

for(i in c(1:10)){
    vioplot2(log_thetas[[i]],
             at = ifelse(i %% 2 == 1, (i + 1) / 2, i / 2),
             side = ifelse(i %% 2 == 1, "left", "right"),
             col = ifelse(i %% 2 == 1, asex_blue, sex_red),
             add = T, h = 0.2)
}
# sex_legend(cex = 1.3)
dev.off()

quartz(type = 'pdf', file = 'figures/tim_w100000_theta_violin_zoomed.pdf')
plot(x=NULL, y=NULL,
     xlim = c(0.5, 5.5), ylim=c(-4, max(unlist(log_thetas))),
     type="n", ann=FALSE, axes=F, bty="n")
axis(1, at=c(1:5), labels = F)
text(c(1:5), par("usr")[3] - 0.2, srt = 10, pos = 1, xpd = TRUE,
     labels = timema_pairs$labels)
axis(2)
mtext(theta_label, side = 2, line = +2, cex = 1.3)

for(i in 1:10){
    to_plot <- log_thetas[[i]][log_thetas[[i]] > -4]
    vioplot2(to_plot,
             at = ifelse(i %% 2 == 1, (i + 1) / 2, i / 2),
             side = ifelse(i %% 2 == 1, "left", "right"),
             col = ifelse(i %% 2 == 1, asex_blue, sex_red),
             add = T, h = 0.2)
}

# sex_legend(cex = 1.3)
dev.off()
asex_indices <- seq(1,10,by=2)

asex_means_pair_ord <- unlist(lapply(thetas[asex_indices], mean))
asex_medians_pair_ord <- unlist(lapply(thetas[asex_indices], median))

#timema_pairs$JC_divergence - divergence between species estimated from transcriptomes

pdf('figures/heterozygosity_mean_vs_age.pdf')
par(mar=c(5.1,4.2,0,2))
plot(asex_means_pair_ord ~ timema_pairs$JC_divergence,
     pch = 20, xlab = 'Mean species pair divergence [%]', ylab = 'mean theta [2Tmu]',
     cex.lab = 1.3, cex = 1.5, bty = 'n', xlim = c(0.7, 1.95)) #, , ylim = c(15, 21.1)
     text(asex_means_pair_ord ~ timema_pairs$JC_divergence,
         labels = timemas$labels[asex_indices], pos = 3)
dev.off()

pdf('figures/heterozygosity_median_vs_age.pdf')
par(mar=c(5.1,4.2,0,2))
plot(asex_medians_pair_ord ~ timema_pairs$JC_divergence,
     pch = 20, xlab = 'Mean species pair divergence [%]', ylab = 'median theta [2Tmu]',
     cex.lab = 1.3, cex = 1.5, bty = 'n', xlim = c(0.7, 1.95)) #, ylim = c(15, 21.1)
     text(asex_medians_pair_ord ~ timema_pairs$JC_divergence,
          labels = timemas$labels[asex_indices], pos = 3)
dev.off()

# #2.8 × 10(-9) / site / generation
# 2.8 * 10e-3 / site / milion years
# 0.00028 - 0.0052
