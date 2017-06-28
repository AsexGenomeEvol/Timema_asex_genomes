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

#1_Tdi 2_Tcm 2_Tsi 3_Tms 4_Tte 5_Tge;
for(i in c(1,3,4,5,7,9)){
    # non-filtered all data
    atlasData[[i]] <- load_thetas(i)
    # filtered theta estimates
    logThetas[[i]] <- log10(filter_theta(atlasData[[i]])[,'theta_MLE'])
}

### violin plot section
    # logThetas 1 3 5 7 9 - asex, 2 4 6 8 10 - sex

require(vioplot)
require(digest)
library(sm)
source('scripts/external/vioplot2.R')
theta_label <- expression(paste("Heterozygosity of 1kbp windows [" , log[10], " " , theta, ']'))

pdf('variant_stats/figures/theta_log_violin.pdf')
plot(x=NULL, y=NULL,
     xlim = c(0.5, 5.5), ylim=c(min(unlist(logThetas)), max(unlist(logThetas))),
     type="n", ann=FALSE, axes=F, bty="n")
axis(1, at=c(1:5), labels = F)
text(c(1:5), par("usr")[3] - 0.2, srt = 10, pos = 1, xpd = TRUE,
     labels = c(expression(italic("T. douglasi / T. poppensis")),
                expression(italic("T. shepardi / T. californicum")),
                expression(italic("T. monikensis / T. cristinae")),
                expression(italic("T. tahoe / T. barmani")),
                expression(italic("T. genevieve / T. podura"))))
axis(2)
mtext(theta_label, side = 2, line = +2, cex = 1.3)

for(i in c(1,3,4,5,7,9)){
    vioplot2(logThetas[[i]],
             at = ifelse(i %% 2 == 1, (i + 1) / 2, i / 2),
             side = ifelse(i %% 2 == 1, "left", "right"),
             col = ifelse(i %% 2 == 1, asex_blue, sex_red),
             add = T, h = 0.2)
}

# sex_legend(cex = 1.3)
dev.off()

pdf('variant_stats/figures/theta_log_violin_zoomed.pdf')
plot(x=NULL, y=NULL,
     xlim = c(0.5, 5.5), ylim=c(-4, max(unlist(logThetas))),
     type="n", ann=FALSE, axes=F, bty="n")
axis(1, at=c(1:5), labels = F)
text(c(1:5), par("usr")[3] - 0.2, srt = 10, pos = 1, xpd = TRUE,
     labels = c(expression(italic("T. douglasi / T. poppensis")),
                expression(italic("T. shepardi / T. californicum")),
                expression(italic("T. monikensis / T. cristinae")),
                expression(italic("T. tahoe / T. barmani")),
                expression(italic("T. genevieve / T. podura"))))
axis(2)
mtext(theta_label, side = 2, line = +2, cex = 1.3)

for(i in c(1,3,4,5,7,9)){
    to_plot <- logThetas[[i]][logThetas[[i]] > -4]
    vioplot2(to_plot,
             at = ifelse(i %% 2 == 1, (i + 1) / 2, i / 2),
             side = ifelse(i %% 2 == 1, "left", "right"),
             col = ifelse(i %% 2 == 1, asex_blue, sex_red),
             add = T, h = 0.2)
}

# sex_legend(cex = 1.3)
dev.off()

human_theta <- read.table('data/pull_data/Theta_estimatesWC1_EM_1kbp_theta_estimates.txt', header = T)
human_yellow <- "#F9EE567D"
human_theta <- filter_theta(human_theta)

plot_Tcm <- function(){
    hist(logThetas[[4]], breaks = 40, col = sex_red, freq = F, xaxt="n",
         xlab = theta_label, main = NA, ylim = c(0, 2.4), cex.lab = 1.3)
    at.x <- outer(1, 10^(c(-8, -3, -2, -1)))
    lab.x <- log10(at.x)
    axis(1, at=lab.x, labels=at.x, las=1, cex.lab = 1.3)
}


plot_Tsi <- function(){
    hist(logThetas[[3]], breaks = 40, col = asex_blue, add = T, freq = F)
}

par(mar=c(6.1,4.2,0,0))
#sets the bottom, left, top and right, def (mar=c(5.1,4.1,4.1,2.1)

quartz(type = 'pdf', file = 'variant_stats/figures/2_Tcm_theta_log_hist.pdf')
# png('variant_stats/figures/2_Tcm_theta_log_hist.png')
    par(mar=c(5.1,4.2,0,2))
    plot_Tcm()
    legend('topright', pch = 20, bty = 'n',
           col = c(sex_red), cex = 1.4,
           expression(italic("T. californicum"~"\u2640"~"\u2642")))
dev.off()


quartz(type = 'pdf', file = 'variant_stats/figures/2_Tcm_Tsi_theta_log_hist.pdf')
    par(mar=c(5.1,4.2,0,2))
    plot_Tcm()
    plot_Tsi()
    legend('topright', pch = 20, bty = 'n',
           col = c(sex_red, asex_blue), cex = 1.4,
           c(expression(italic("T. californicum"~"\u2640"~"\u2642")),
             expression(italic("T. shepardi"~"\u2640"))))
dev.off()

quartz(type = 'pdf', file = 'variant_stats/figures/2_Tcm_Tsi_theta_log_hist_with_human.pdf')
    par(mar=c(5.1,4.2,0,2))
    plot_Tcm()
    plot_Tsi()
    hist(log10(human_theta$theta_MLE), breaks = 40, col = human_yellow, add = T, freq = F)
    legend('topright', pch = 20, bty = 'n', cex = 1.4,
           col = c(sex_red, asex_blue, human_yellow),
           c(expression(italic("T. californicum"~"\u2640"~"\u2642")),
             expression(italic("T. shepardi"~"\u2640")),
             expression(italic("H. sapiens"~"\u2640"~"\u2642"))))
dev.off()

means <- list()
for(i in c(1,3,4,5,7,9)){
    means[[i]] <- mean(10^(logThetas[[i]]))
}

#plot(unlist(means[seq(1,10,by=2)]) ~ divergence_times, pch = 20)
# unlist(lapply(logThetas[seq(1,10,by=2)], function(x){ sum(x < -5)})) / 1000
pdf('variant_stats/figures/heterozygosity_vs_age.pdf')
par(mar=c(5.1,4.2,0,2))
heterozygous_portion <- unlist(lapply(logThetas[seq(1,10,by=2)], function(x){ mean(x > -5)})) * 100
plot(heterozygous_portion ~ divergence_times,
     pch = 20, xlab = 'Mean species pair divergence [%]', ylab = 'proportion of heterozygous regions in genome [%]',
     cex.lab = 1.3, cex = 1.5, bty = 'n', xlim = c(0.7, 1.95), ylim = c(15, 21.1))
text(heterozygous_portion ~ divergence_times,
    labels = timema_labels[c(1,3,5,7,9)], pos = 3)
dev.off()

pdf('variant_stats/figures/mean_heterozygosity_vs_age.pdf')
plot(unlist(lapply(logThetas[seq(1,10,by=2)], function(x){ mean(10^x)})) ~ divergence_times,
     pch = 20, xlab = 'Mean species pair divergence [%]', ylab = expression(paste('Mean heterozygosity [', theta, ']')),
     cex.lab = 1.3, cex = 1.5, bty = 'n', xlim = c(0.7,1.95), ylim = c(0.00095, 0.00120))
text(unlist(lapply(logThetas[seq(1,10,by=2)], function(x){ mean(10^x)})) ~ divergence_times,
    labels = timema_labels[c(1,3,5,7,9)], pos = 3)
dev.off()

pdf('variant_stats/figures/mean_heterozygosity_vs_very_high_estimates.pdf')
plot(unlist(lapply(logThetas[seq(1,10,by=2)], function(x){ mean(10^x)})) ~ unlist(lapply(logThetas, function(x){ sum(x > -1.5)})[seq(1,10,by = 2)]),
     pch = 20, xlab = '# estimates > 0.03', ylab = expression(paste('Mean heterozygosity [', theta, ']')),
     cex.lab = 1.3, cex = 1.5, bty = 'n', ylim = c(0.00095, 0.00120))
dev.off()
# > 0.3



# #2.8 × 10(-9) / site / generation
# 300*10^6
# 2.8 * 1.3 * 2 * 10e6
# # 7 280 000
# # 0.00028 after 0.5 * 10e6
# # 0.00056 after 1 * 10e6
# # 0.00112 after 2 * 10e6
# 0.00095 - 0.00120

##### COVERAGE PLOTS
i <- 4
head(atlasData[[i]])

sex_coverages <- atlasData[[i]]$coverage[atlasData[[i]]$coverage < 100]
sex_mid <- median(sex_coverages)
sex_min <- sex_mid * 0.5
sex_max <- sex_mid * 1.5

asex_coverages <- atlasData[[i-1]]$coverage[atlasData[[i-1]]$coverage < 100]
asex_mid <- median(asex_coverages)
asex_min <- asex_mid * 0.5    #nrow(sp_data)
asex_max <- asex_mid * 1.5

quartz(type = 'pdf', file = 'variant_stats/figures/2_Tcm_Tsi_coverage_hist.pdf')
    par(mar=c(5.1,4.2,2,2))
hist(sex_coverages, xlim = c(0, 100), breaks = 50, col = sex_red, ylim = c(0,340000),
    xlab = 'Coverge', main = '', cex.lab = 1.3)
#lines(c(sex_mid, 0), c(sex_mid, 10e6), col = sex_red, lty = 1)
hist(asex_coverages, xlim = c(0, 100), breaks = 50, col = asex_blue, add = T)
#lines(c(asex_mid, 0), c(asex_mid, 10e6), col = asex_blue, lty = 1)
legend('topright', pch = 20, bty = 'n',
       col = c(sex_red, asex_blue), cex = 1.4,
       c(expression(italic("T. californicum"~"\u2640"~"\u2642")),
         expression(italic("T. shepardi"~"\u2640"))))
# lines(c(sex_min, 0), c(sex_min, 10e6), col = sex_red, lty = 1, lwd = 2)
# lines(c(sex_max, 0), c(sex_max, 10e6), col = sex_red, lty = 1, lwd = 2)
# lines(c(asex_min, 0), c(asex_min, 10e6), col = asex_blue, lty = 1, lwd = 2)
# lines(c(asex_max, 0), c(asex_max, 10e6), col = asex_blue, lty = 1, lwd = 2)
dev.off()

# # only reasonable coverages
# hist(Tcm$coverage)
# # no correlation between coverage and theta estimate
# plot(Tcm$theta_MLE ~ Tcm$coverage, pch = 20)
# # likelihoods seems to be of one distribution
# hist(Tcm$LL)
#hist(Tcm$theta_C95_u - Tcm$theta_C95_l, breaks = 2000, xlim = c(0,0.1))
