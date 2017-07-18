# load_thetas.R from AsexStats
# plot_log_hist.R
# sex_legend.R

Os_thetas <- load_thetas('1_Os', 'ref', '1_Os_b1v01', 10000)
Os_thetas <- Os_thetas[(Os_thetas$end - Os_thetas$start) > 5000,]
Os_thetas <- Os_thetas[!is.na(Os_thetas$theta_MLE),]

On_thetas <- load_thetas('1_On', 'ref', '1_On_b1v01', 10000)
On_thetas <- On_thetas[(On_thetas$end - On_thetas$start) > 5000,]
On_thetas <- On_thetas[!is.na(On_thetas$theta_MLE),]

Os_logthetas <- log10(Os_thetas$theta_MLE)
On_logthetas <- log10(On_thetas$theta_MLE)

theta_log_label <- expression(paste("Heterozygosity of 10kbp windows [" , log[10], " " , theta, ']'))

quartz(type = 'pdf', file = 'figures/Os_On_theta_loghist.pdf')
hist(Os_logthetas, breaks = 40, col = sex_red, freq = F, xaxt="n",
     xlab = theta_log_label, ylim = c(0, 2.4), cex.lab = 1.3,
     main = 'Os / On')
at.x <- outer(1, 10^(c(-8, -3, -2, -1)))
lab.x <- log10(at.x)
axis(1, at=lab.x, labels=at.x, las=1, cex.lab = 1.3)
hist(On_logthetas, breaks = 40, col = asex_blue, add = T, freq = F)
sex_legend()
dev.off()

quartz(type = 'pdf', file = 'figures/Os_On_theta_ylog_hist.pdf')
plot_log_hist(sex_data = Os_thetas$theta_MLE, asex_data = On_thetas$theta_MLE,
              breaks = 160, barwidth = 4,
              xlab = expression(paste("Heterozygosity of 10kbp windows [", theta, ']')))
sex_legend()
dev.off()

summary(Os_thetas$theta_MLE)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.000000 0.001410 0.004713 0.005692 0.008432 0.155200

summary(On_thetas$theta_MLE)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.000000 0.004383 0.007777 0.007986 0.010590 0.115100

# Oppiella nova / On6 (asexual)
# Oppiella subpectinata / Os5 (sexual)
# Atropacarus striculus / As6 (asexual)
# Steganacarus magnus / Sm2 (sexual)



file <- 'data/2_Sm/variant_calls/ref/atlas/ref_to_2_Sm_b1v01_w10000_theta_estimates_killed.txt'
Sm_thetas <- read.table(file, header = T)
Sm_thetas <- Sm_thetas[(Sm_thetas$end - Sm_thetas$start) > 5000,]
Sm_thetas <- Sm_thetas[!is.na(Sm_thetas$theta_MLE),]

file <- 'data/2_As/variant_calls/ref/atlas/ref_to_2_As_b1v01_w10000_theta_estimates.txt'
As_thetas <-  read.table(file, header = T)
As_thetas <- As_thetas[(As_thetas$end - As_thetas$start) > 5000,]
As_thetas <- As_thetas[!is.na(As_thetas$theta_MLE),]

Sm_logthetas <- log10(Sm_thetas$theta_MLE)
As_logthetas <- log10(As_thetas$theta_MLE)

theta_log_label <- expression(paste("Heterozygosity of 10kbp windows [" , log[10], " " , theta, ']'))

quartz(type = 'pdf', file = 'figures/Sm_As_theta_loghist.pdf')
hist(Sm_logthetas, breaks = 40, col = sex_red, freq = F, xaxt="n",
     xlab = theta_log_label, ylim = c(0, 2.4), cex.lab = 1.3,
     main = 'Sm / As')
at.x <- outer(1, 10^(c(-8, -3, -2, -1)))
lab.x <- log10(at.x)
axis(1, at=lab.x, labels=at.x, las=1, cex.lab = 1.3)
hist(As_logthetas, breaks = 40, col = asex_blue, add = T, freq = F)
sex_legend()
dev.off()

quartz(type = 'pdf', file = 'figures/Sm_As_theta_ylog_hist.pdf')
plot_log_hist(sex_data = Sm_thetas$theta_MLE, asex_data = As_thetas$theta_MLE,
              breaks = 160, barwidth = 4,
              xlab = expression(paste("Heterozygosity of 10kbp windows [", theta, ']')))
sex_legend()
dev.off()

#### SEX SEX figures

dark_sex = "#C95048ED"
light_sex = "#E9A070B5"
plot_log_hist(sex_data = Sm_thetas$theta_MLE, asex_data = On_thetas$theta_MLE,
              breaks = 160, barwidth = 4, col = c(dark_sex, light_sex),
              xlab = expression(paste("Heterozygosity of 10kbp windows [", theta, ']')))

# Polite mites

##### READ_DATA
args <- commandArgs(trailingOnly=TRUE)

library(AsexStats)

ind <- args[1]
ref <- args[2]
window <- args[3]
ind <- 'ref'
ref <- 'b1v01'
window <- '10000'

thetas <- list()
logThetas <- list()
atlasData <- list()

mites <- list()
              #   ASEX    SEX
mites$codes <- c("1_On","1_Os",
                 "2_As","2_Sm")

files <- c(
    'data/1_On/variant_calls/ref/atlas/ref_to_1_On_b1v01_w10000_theta_estimates.txt',
    'data/1_Os/variant_calls/ref/atlas/ref_to_1_Os_b1v01_w10000_theta_estimates.txt',
    'data/2_As/variant_calls/ref/atlas/ref_to_2_As_b1v01_w10000_theta_estimates.txt',
    'data/2_Sm/variant_calls/ref/atlas/ref_to_2_Sm_b1v01_w10000_theta_estimates_killed.txt')

for(i in 1:4){
    # non-filtered all data
    atlasData[[i]] <- read.table(files[i], header = T)
    # filtered theta estimates
    thetas[[i]] <- filter_thetas(atlasData[[i]], filt_cov = F, window_size = 9999)[,'theta_MLE']
    logThetas[[i]] <- log10(thetas[[i]])
}

library(sm)
source('/Volumes/dump/projects/timema/AsexStats/R/vioplot2.R')

theta_label <- expression(paste("Heterozygosity of 10kbp windows [" , log[10], " " , theta, ']'))

quartz(type = 'pdf', file = 'figures/mites_w10000_theta_violin.pdf')
plot(x=NULL, y=NULL,
     xlim = c(0.5, 2.5), ylim=c(min(unlist(logThetas)), max(unlist(logThetas))),
     type="n", ann=FALSE, axes=F, bty="n")
axis(1, at=c(1:5), labels = F)
text(c(1:5), par("usr")[3] - 0.2, srt = 10, pos = 1, xpd = TRUE,
     labels = c('On / Os', 'As / Sm'))
axis(2)
mtext(theta_label, side = 2, line = +2, cex = 1.3)

for(i in 1:4){
    vioplot2(logThetas[[i]],
             at = ifelse(i %% 2 == 1, (i + 1) / 2, i / 2),
             side = ifelse(i %% 2 == 1, "left", "right"),
             col = ifelse(i %% 2 == 1, asex_blue, sex_red),
             add = T, h = 0.2)
}
# sex_legend(cex = 1.3)
dev.off()
