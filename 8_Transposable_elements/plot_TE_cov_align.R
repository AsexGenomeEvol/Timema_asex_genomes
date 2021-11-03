library(plyr)
library(stringr)
library(ggplot2)
library(cowplot)


### read data (large files! take a while to read)
#setwd("/Users/dparker/Documents/University/Lausanne/Timema_genomes/TEs/")

dat1_Tdi <-  read.csv("1_Tdi_TE_cov_align.txt", sep = "\t")
dat1_Tps <-  read.csv("1_Tps_TE_cov_align.txt", sep = "\t")
dat1_Tcm <-  read.csv("2_Tcm_TE_cov_align.txt", sep = "\t")
dat1_Tsi <-  read.csv("2_Tsi_TE_cov_align.txt", sep = "\t")
dat1_Tce <-  read.csv("3_Tce_TE_cov_align.txt", sep = "\t")
dat1_Tms <-  read.csv("3_Tms_TE_cov_align.txt", sep = "\t")
dat1_Tbi <-  read.csv("4_Tbi_TE_cov_align.txt", sep = "\t")
dat1_Tte <-  read.csv("4_Tte_TE_cov_align.txt", sep = "\t")
dat1_Tge <-  read.csv("5_Tge_TE_cov_align.txt", sep = "\t")
dat1_Tpa <-  read.csv("5_Tpa_TE_cov_align.txt", sep = "\t")

## plotting function

plot_TE_cov <- function(df) {
  
  title_gg = paste(unique(df$sp)[1], unique(df$sp)[2])
  print( title_gg )
  df$rep_mode <-   ifelse(df$sp == "Tps", "sex", 
                          ifelse(df$sp == "Tdi", "asex",
                                 ifelse(df$sp == "Tcm", "sex",
                                        ifelse(df$sp == "Tsi", "asex",
                                               ifelse(df$sp == "Tce", "sex",
                                                      ifelse(df$sp == "Tms", "asex",
                                                             ifelse(df$sp == "Tbi", "sex",
                                                                    ifelse(df$sp == "Tte", "asex",
                                                                           ifelse(df$sp == "Tpa", "sex",
                                                                                  ifelse(df$sp == "Tge", "asex", "ERROR"))))))))))

p1 <- ggplot(df, aes(x=div_a, r_cov_n, col = rep_mode)) +
  theme_bw() +
  geom_line() +
  xlab("Kimura substitution level (CpG adjusted)") + 
  ylab("Normalised read coverage") + 
  scale_color_manual(values = c("asex" = "blue", "sex" = "red")) + ggtitle(title_gg)

return(p1)
}


####### specify which TEs I want to keep (I.e. all but unknown and chimeric TEs)

TE_keep_list <- c(
'DNA/Academ', 'DNA/CMC', 'DNA/Crypton', 'DNA/Mite', 'DNA/Chapaev', 'DNA/Ginger',
'DNA/Harbinger', 'DNA/hAT', 'DNA/Kolobok', 'DNA/Maverick', 'DNA',  'DNA/Merlin', 
'DNA/MULE', 'DNA/P', 'DNA/PiggyBac', 'DNA/Sola', 'DNA/TcMar', 'DNA/Transib', 
'DNA/Zator', 'DNA/Dada', 'DNA/MuDr', 'DNA/MITE', 'DNA/ISL2EU', 'DIRS/Other',
'RC/Helitron', 'LTR/BEL', 'LTR/DIRS', 'LTR/Ngaro', 'LTR/Pao', 'LTR/Copia',
'LTR/Gypsy', 'LTR/ERVL','LTR', 'LTR/ERV1', 'LTR/ERV', 'LTR/ERVK', 
'LTR/Penelope', 'LINE/L1', 'LINE', 'LINE/RTE', 'LINE/CR1', 'LINE/Rex-Babar','LINE/L2', 
'LINE/Proto2', 'LINE/LOA', 'LINE/R1', 'LINE/I', 'LINE/Jockey', 'LINE/Jockey-I',
'LINE/Dong-R4', 'LINE/R2', 'LINE/Penelope', 'LINE/Poseidon', 'LINE/CRE', 'LINE/Tx1',
'LINE/Crack', 'Retroposon/SVA', 'SINE', 'SINE/5S', 'SINE/7SL', 'SINE/Alu', 
'SINE/tRNA', 'SINE/tRNA-Alu', 'SINE/tRNA-RTE', 'SINE/RTE', 'SINE/Deu',
'SINE/tRNA-V', 'SINE/MIR', 'SINE/U', 'SINE/tRNA-7SL', 'SINE/tRNA-CR1')


calc_cov_by_div_TEs_together_2 <- function(df, sp){
  ## make bins
  df$div_g <- data.frame(df$Kimura_with_divCpGMod,bin=cut(df$Kimura_with_divCpGMod,seq(0,100, 2),include.lowest=TRUE))$bin
  df2 <- na.omit(df)
  
  ### drop unknown and chimeric TEs # note comment out these lines if want to include everything. It produces similar results.
  
  df2_excl <- subset(df2, ! TE_A_class_s %in% TE_keep_list)
  df2      <- subset(df2,   TE_A_class_s %in% TE_keep_list)
  
  ## print TEs used/discarded as a check
  print("TEs kept:")
  print(levels(as.factor(df2$TE_A_class_s)))  
  print("TEs excluded:")
  print(levels(as.factor(df2_excl$TE_A_class_s)))  

  ## all TEs together
  df3 <- ddply(df2, .(div_g), summarize,  Nreadcounts=sum(Nreadcounts), feat_len=sum(feat_len))
  df3$r_cov <- df3$Nreadcounts / df3$feat_len * 100
  
  ## better bin name
  df3$div_a <- as.numeric(str_split_fixed(gsub("\\[", "", gsub("\\(", "", as.character(df3$div_g))), ",", 2)[,1]) + 1
  
  ### set max div 35 as few sites over this div
  df3 <- subset(df3, df3$div_a <= 35)
  
  ### adjust for overall differences in coverage
  df3$r_cov_n <- df3$r_cov / median(df3$r_cov)
  
  ## add sp
  df3$sp <- rep(sp, length(df3[,1]))
  
  return(df3)  
}

Tdi_TE_cov_2 <- calc_cov_by_div_TEs_together_2(dat1_Tdi, "Tdi")
Tps_TE_cov_2 <- calc_cov_by_div_TEs_together_2(dat1_Tps, "Tps")

Tcm_TE_cov_2 <- calc_cov_by_div_TEs_together_2(dat1_Tcm, "Tcm")
Tsi_TE_cov_2 <- calc_cov_by_div_TEs_together_2(dat1_Tsi, "Tsi")

Tce_TE_cov_2 <- calc_cov_by_div_TEs_together_2(dat1_Tce, "Tce")
Tms_TE_cov_2 <- calc_cov_by_div_TEs_together_2(dat1_Tms, "Tms")

Tbi_TE_cov_2 <- calc_cov_by_div_TEs_together_2(dat1_Tbi, "Tbi")
Tte_TE_cov_2 <- calc_cov_by_div_TEs_together_2(dat1_Tte, "Tte")

Tpa_TE_cov_2 <- calc_cov_by_div_TEs_together_2(dat1_Tpa, "Tpa")
Tge_TE_cov_2 <- calc_cov_by_div_TEs_together_2(dat1_Tge, "Tge")


TpsTdi_TE_cov_2 <- rbind(Tps_TE_cov_2, Tdi_TE_cov_2)
TcmTsi_TE_cov_2 <- rbind(Tcm_TE_cov_2, Tsi_TE_cov_2)
TceTms_TE_cov_2 <- rbind(Tce_TE_cov_2, Tms_TE_cov_2)
TbiTte_TE_cov_2 <- rbind(Tbi_TE_cov_2, Tte_TE_cov_2)
TpaTge_TE_cov_2 <- rbind(Tpa_TE_cov_2, Tge_TE_cov_2)

## I provide these files in for convenience in 8_Transposable_elements/R/Repeatmasker/coverage
write.csv(TpsTdi_TE_cov_2, "TpsTdi_TE_cov_align.csv")
write.csv(TcmTsi_TE_cov_2, "TcmTsi_TE_cov_align.csv")
write.csv(TceTms_TE_cov_2, "TceTms_TE_cov_align.csv")
write.csv(TbiTte_TE_cov_2, "TbiTte_TE_cov_align.csv")
write.csv(TpaTge_TE_cov_2, "TpaTge_TE_cov_align.csv")

# ### read back in 
# TpsTdi_TE_cov_2 <- read.csv("TpsTdi_TE_cov_align.csv")
# TcmTsi_TE_cov_2 <- read.csv("TcmTsi_TE_cov_align.csv")
# TceTms_TE_cov_2 <- read.csv("TceTms_TE_cov_align.csv")
# TbiTte_TE_cov_2 <- read.csv("TbiTte_TE_cov_align.csv")
# TpaTge_TE_cov_2 <- read.csv("TpaTge_TE_cov_align.csv")
# 


### plot as line plot

pdf("TE_cov_align_lineplot.pdf", width = 7, height = 8)
plot_grid(
  plot_TE_cov(TpsTdi_TE_cov_2 ),
  plot_TE_cov(TcmTsi_TE_cov_2 ),
  plot_TE_cov(TceTms_TE_cov_2 ),
  plot_TE_cov(TbiTte_TE_cov_2),
  plot_TE_cov(TpaTge_TE_cov_2), ncol = 2)
dev.off()
getwd() ## where has my plot gone....?

png(filename = "TE_cov_align_lineplot.png", width = 7, height = 8, units = "in", bg = "white", res = 300)
plot_grid(
  plot_TE_cov(TpsTdi_TE_cov_2 ),
  plot_TE_cov(TcmTsi_TE_cov_2 ),
  plot_TE_cov(TceTms_TE_cov_2 ),
  plot_TE_cov(TbiTte_TE_cov_2),
  plot_TE_cov(TpaTge_TE_cov_2), ncol = 2)
dev.off()
getwd() ## where has my plot gone....



### plot as bar plot

plot_TE_cov_bar <- function(df) {
  
  ## better bin name
  df$div_a1   <- paste(df$div_a - 1, "-",   df$div_a + 1, sep = "")
  df$div_a1_o <- ordered(as.factor(df$div_a1), levels = c('0-2', '2-4', '4-6', '6-8', '8-10', '10-12', '12-14', '14-16', '16-18', '18-20', '20-22', '22-24', '24-26', '26-28', '28-30', '30-32', '32-34', '34-36'))
  print(str(df))
  print(df)
  
  title_gg = paste(unique(df$sp)[1], unique(df$sp)[2])
  print( title_gg )
  df$rep_mode <-   ifelse(df$sp == "Tps", "sex", 
                          ifelse(df$sp == "Tdi", "asex",
                                 ifelse(df$sp == "Tcm", "sex",
                                        ifelse(df$sp == "Tsi", "asex",
                                               ifelse(df$sp == "Tce", "sex",
                                                      ifelse(df$sp == "Tms", "asex",
                                                             ifelse(df$sp == "Tbi", "sex",
                                                                    ifelse(df$sp == "Tte", "asex",
                                                                           ifelse(df$sp == "Tpa", "sex",
                                                                                  ifelse(df$sp == "Tge", "asex", "ERROR"))))))))))
  
  p1 <- ggplot(df, aes(x=div_a1_o, r_cov_n, fill = rep_mode)) +
    theme_bw() +
    geom_bar(stat="identity", position="dodge") +
    xlab("Kimura substitution level (CpG adjusted)") + 
    ylab("Normalised read coverage") + 
    scale_fill_manual(values = c("asex" = "#92C5DECD", "sex" =  "#D6604DED")) + ggtitle(title_gg)
  
  p2 <- p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(p2)
}


pdf("TE_cov_align_barplot.pdf", width = 9, height = 8)
plot_grid(
  plot_TE_cov_bar(TpsTdi_TE_cov_2 ),
  plot_TE_cov_bar(TcmTsi_TE_cov_2 ),
  plot_TE_cov_bar(TceTms_TE_cov_2 ),
  plot_TE_cov_bar(TbiTte_TE_cov_2),
  plot_TE_cov_bar(TpaTge_TE_cov_2), ncol = 2)
dev.off()
getwd() ## where has my plot gone....?

png(filename = "TE_cov_align_barplot.png", width = 9, height = 8, units = "in", bg = "white", res = 300)
plot_grid(
  plot_TE_cov_bar(TpsTdi_TE_cov_2 ),
  plot_TE_cov_bar(TcmTsi_TE_cov_2 ),
  plot_TE_cov_bar(TceTms_TE_cov_2 ),
  plot_TE_cov_bar(TbiTte_TE_cov_2),
  plot_TE_cov_bar(TpaTge_TE_cov_2), ncol = 2)
dev.off()
getwd() ## where has my plot gone....






