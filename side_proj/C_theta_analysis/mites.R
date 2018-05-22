library(AsexStats)

mites <- list()
              #   ASEX    SEX
mites$labels <- c("1_On","1_Os",
                  "2_As","2_Sm")

ref_files <- c(
    'data/1_On/variant_calls/ref/atlas/ref_to_1_On_b1v01_w10000_theta_estimates.txt',
    'data/1_Os/variant_calls/ref/atlas/ref_to_1_Os_b1v01_w10000_theta_estimates.txt',
    'data/2_As/variant_calls/ref/atlas/ref_to_2_As_b1v01_w10000_theta_estimates.txt',
    'data/2_Sm/variant_calls/ref/atlas/ref_to_2_Sm_b1v01_w10000_theta_estimates_killed.txt')
reseq_files <- sapply(mites$labels, function(x) { dir(paste0("data/", x), pattern = "theta_estimates.txt", full.names = T) })
files <- as.vector(rbind(ref_files, reseq_files))
atlas_data <- lapply(files, read.table, header = T)
thetas <- lapply(atlas_data, function(x){filter_thetas(x, filt_cov = F, window_size = 9999)[,'theta_MLE']} )
log_thetas <- lapply(thetas, log10)

theta_log_label <- expression(paste("Heterozygosity of 10kbp windows [" , log[10], " " , theta, ']'))

### On / Os pair
## x log hist
quartz(type = 'pdf', file = 'figures/Os_On_theta_loghist.pdf')
hist(log_thetas[[2]], breaks = 40, col = sex_red, freq = F, xaxt="n",
     xlab = theta_log_label, ylim = c(0, 2.4), cex.lab = 1.3,
     main = 'On / Os')
at.x <- outer(1, 10^(c(-8, -3, -2, -1)))
lab.x <- log10(at.x)
axis(1, at=lab.x, labels=at.x, las=1, cex.lab = 1.3)
hist(log_thetas[[1]], breaks = 40, col = asex_blue, add = T, freq = F)
sex_legend()
dev.off()

## y log hist
quartz(type = 'pdf', file = 'figures/Os_On_theta_ylog_hist.pdf')
plot_log_hist(sex_data = thetas[[2]], asex_data = thetas[[1]],
              breaks = 160, barwidth = 3, cex.axis = 1.3,
              xlab = expression(paste("Heterozygosity of 10kbp windows [", theta, ']')))
sex_legend()
dev.off()

### As / Sm pair
## x log hist
quartz(type = 'pdf', file = 'figures/Sm_As_theta_loghist.pdf')
hist(log_thetas[[4]], breaks = 40, col = sex_red, freq = F, xaxt="n",
     xlab = theta_log_label, ylim = c(0, 2.4), cex.lab = 1.3,
     main = 'As / Sm')
at.x <- outer(1, 10^(c(-8, -3, -2, -1)))
lab.x <- log10(at.x)
axis(1, at=lab.x, labels=at.x, las=1, cex.lab = 1.3)
hist(log_thetas[[3]], breaks = 40, col = asex_blue, add = T, freq = F)
sex_legend()
dev.off()

## y log hist
quartz(type = 'pdf', file = 'figures/Sm_As_theta_ylog_hist.pdf')
plot_log_hist(asex_data = thetas[[3]], sex_data = thetas[[4]],
              breaks = 40, barwidth = 16, cex.axis = 1.3,
              xlab = expression(paste("Heterozygosity of 10kbp windows [", theta, ']')))
sex_legend()
dev.off()

## summary table
theta_summaries <- lapply(thetas, summary)
get_stat <- function(input_list, stat_name, round_to = 4){
    round(unlist(lapply(input_list, function(x){ x[stat_name] })), round_to)
}

mite_summary <- data.frame(mite = substr(files, 11, 13),
                           mean = get_stat(theta_summaries, 'Mean'),
                           st_quantile = get_stat(theta_summaries, '1st Qu.'),
                           median = get_stat(theta_summaries, 'Median'),
                           rd_quantile = get_stat(theta_summaries, '3rd Qu.'))
write.table(mite_summary, 'stats/ref_theta_10k_summary.tsv',
            quote = F, sep = '\t', row.names = F)

quartz(type = 'pdf', file = 'figures/mites_heterozygosity_mediabs.pdf')
    bp_data <- matrix(mite_summary$median, ncol = 4)
    bp_lq <- matrix(mite_summary$st_quantile, ncol = 4)
    bp_uq <- matrix(mite_summary$rd_quantile, ncol = 4)
    bp <- barplot(bp_data, names = mites$labels, ylab = 'median Heterozygosity',
                  col = rep(c(asex_blue, sex_red), each = 6), beside = TRUE, ylim = c(0,1))

    segments(bp, bp_data - bp_lq, bp, bp_data + bp_uq, lwd = 1.5)
    arrows(bp, bp_data - bp_lq, bp, bp_data + bp_uq, lwd = 1.5, angle = 90,
           code = 3, length = 0.05)
dev.off()

#### SEX SEX figures
# dark_sex = "#C95048ED"
# light_sex = "#E9A070B5"
# plot_log_hist(sex_data = Sm_thetas$theta_MLE, asex_data = On_thetas$theta_MLE,
#               breaks = 160, barwidth = 4, col = c(dark_sex, light_sex),
#               xlab = expression(paste("Heterozygosity of 10kbp windows [", theta, ']')))

## Violiolin plots

log_thetas <- lapply(log_thetas, function(x) x[!is.infinite(x)])
log_thetas <- lapply(log_thetas, function(x) x[!is.na(x)])

quartz(type = 'pdf', file = 'figures/all_mites_w10000_theta_violin.pdf')
plot(x=NULL, y=NULL,
     xlim = c(0.5, 2.5), ylim = range(log_thetas, na.rm = T, finite = T),
     type="n", ann=FALSE, axes=F, bty="n")
axis(1, at=c(1:5), labels = F)
text(c(1:5), par("usr")[3] - 0.2, srt = 0, pos = 1, xpd = TRUE,
     labels = c('On / Os', 'As / Sm'))
axis(2)
mtext(theta_log_label, side = 2, line = +2, cex = 1.3)

for(i in 1:24){
    vioplot2(log_thetas[[i]],
             at = ceiling(1:24 / 12)[i],
             side = ifelse(ceiling(1:24 / 6)[i] %% 2 == 1, "left", "right"),
             col = ifelse(ceiling(1:24 / 6)[i] %% 2 == 1, asex_blue, sex_red),
             add = T, h = 0.2)
}
# sex_legend(cex = 1.3)
dev.off()
