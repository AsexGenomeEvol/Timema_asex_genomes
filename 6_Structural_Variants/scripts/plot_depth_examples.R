library('RColorBrewer')
library('AsexStats')
library('grDevices')

# LOAD SNPs
tab_filenames <- paste0('data/SNP_calls/', c('1_Tdi', '1_Tps') ,'_reduced_variants.tsv')
variant_tab <- lapply(tab_filenames, read.table, stringsAsFactors = F, col.names = c('scf', 'pos', 'qual', paste0('g', 1:5), paste0('d', 1:5)))

asex_pal <- brewer.pal(3, 'Blues')[c(3, 2)]
asex_pal[2] <- adjustcolor(asex_pal[2], alpha.f = 0.6)
sex_pal <- brewer.pal(3, 'Reds')[c(3, 2)]
sex_pal[2] <- adjustcolor(sex_pal[2], alpha.f = 0.6)

Tdi_01 <- data.frame(genotype = variant_tab[[1]]$g1, depth = variant_tab[[1]]$d1)
Tdi_01 <- Tdi_01[Tdi_01$genotype != './.', ]
Tdi_01$depth <- as.numeric(Tdi_01$depth)

Tps_01 <- data.frame(genotype = variant_tab[[2]]$g1, depth = variant_tab[[2]]$d1)
Tps_01 <- Tps_01[Tps_01$genotype != './.', ]
Tps_01$depth <- as.numeric(Tps_01$depth)

### LOAD SVS
source('C_SV_calling/vcf_processing_fctions.R')

sp = '1_Tdi'
sp_short <- substr(sp, 3, 5)
SV_files <- paste0('data/manta_SV_calls/data/', sp, '/', sp_short, '_0', 1 ,'_manta/results/variants/diploidSV_reduced.vcf')

ind = 1

SV_tab <- read.table(SV_files[ind], header = F)
genotypes <- substr(SV_tab$V10, 1, 3)
SR_cov_1 <- sapply(strsplit(SV_tab[,10], ':'), str2depth, 6, 1)
SR_cov_2 <- sapply(strsplit(SV_tab[,10], ':'), str2depth, 6, 2)
SR_cov <- SR_cov_1 + SR_cov_2

SVs_to_plot <- SR_cov < 80
Tdi_SV_covs <- SR_cov[SVs_to_plot]
Tdi_SV_gen <- genotypes[SVs_to_plot]

sp = '1_Tps'
sp_short <- substr(sp, 3, 5)
SV_files <- paste0('data/manta_SV_calls/data/', sp, '/', sp_short, '_0', 1 ,'_manta/results/variants/diploidSV_reduced.vcf')

ind = 1

SV_tab <- read.table(SV_files[ind], header = F)
genotypes <- substr(SV_tab$V10, 1, 3)
SR_cov_1 <- sapply(strsplit(SV_tab[,10], ':'), str2depth, 6, 1)
SR_cov_2 <- sapply(strsplit(SV_tab[,10], ':'), str2depth, 6, 2)
SR_cov <- SR_cov_1 + SR_cov_2

SVs_to_plot <- SR_cov < 80
Tps_SV_covs <- SR_cov[SVs_to_plot]
Tps_SV_gen <- genotypes[SVs_to_plot]

par(mfrow = c(2, 2))

hist(Tdi_01[Tdi_01$genotype == '1/1', 'depth'], col = asex_pal[1], breaks = 60, probability = T, main = 'T. douglasi', xlab = 'coverage')
hist(Tdi_01[Tdi_01$genotype == '0/1', 'depth'], col = asex_pal[2], breaks = 60, add = T, probability = T)
legend('topright', c('1/1', '0/1'), col = asex_pal, bty = 'n', pch = 20)

hist(Tps_01[Tps_01$genotype == '1/1', 'depth'], col = sex_pal[1], breaks = 60, probability = T, main = 'T. poppensis', ylim = c(0, 0.08), xlab = 'coverage' )
hist(Tps_01[Tps_01$genotype == '0/1', 'depth'], col = sex_pal[2], breaks = 120, add = T, probability = T)
legend('topright', c('1/1', '0/1'), col = sex_pal, bty = 'n', pch = 20)

hist(Tdi_SV_covs[Tdi_SV_gen == '1/1'], col = asex_pal[1], breaks = 60, probability = T, xlab = 'coverage', main = '')
hist(Tdi_SV_covs[Tdi_SV_gen == '0/1'], col = asex_pal[2], breaks = 60, add = T, probability = T)
legend('topright', c('1/1', '0/1'), col = asex_pal, bty = 'n', pch = 20)

hist(Tps_SV_covs[Tps_SV_gen == '1/1'], col = sex_pal[1], breaks = 60, probability = T, xlab = 'coverage', main = '')
hist(Tps_SV_covs[Tps_SV_gen == '0/1'], col = sex_pal[2], breaks = 60, add = T, probability = T)
legend('topright', c('1/1', '0/1'), col = sex_pal, bty = 'n', pch = 20)

