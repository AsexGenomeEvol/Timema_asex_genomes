library(ggplot2)
#library(cowplot)
library(ggpubr)

# https://www.r-graph-gallery.com/48-grouped-barplot-with-ggplot2/

setwd("/Users/admin/Documents/schwander/timema/TE/landscapes_knownTEs")

Tps_data <- read.table("1_Tps_b3v07_R_data.txt", sep="\t", header=TRUE, comment.char="#",
                   na.strings=c("","NA"), stringsAsFactors=FALSE,
                   quote="", fill=FALSE)

Tdi_data <- read.table("1_Tdi_b3v07_R_data.txt", sep="\t", header=TRUE, comment.char="#",
                       na.strings=c("","NA"), stringsAsFactors=FALSE,
                       quote="", fill=FALSE)

Tcm_data <- read.table("2_Tcm_b3v07_R_data.txt", sep="\t", header=TRUE, comment.char="#",
                       na.strings=c("","NA"), stringsAsFactors=FALSE,
                       quote="", fill=FALSE)

Tsi_data <- read.table("2_Tsi_b3v07_R_data.txt", sep="\t", header=TRUE, comment.char="#",
                       na.strings=c("","NA"), stringsAsFactors=FALSE,
                       quote="", fill=FALSE)

Tce_data <- read.table("3_Tce_b3v07_R_data.txt", sep="\t", header=TRUE, comment.char="#",
                       na.strings=c("","NA"), stringsAsFactors=FALSE,
                       quote="", fill=FALSE)

Tms_data <- read.table("3_Tms_b3v07_R_data.txt", sep="\t", header=TRUE, comment.char="#",
                       na.strings=c("","NA"), stringsAsFactors=FALSE,
                       quote="", fill=FALSE)

Tbi_data <- read.table("4_Tbi_b3v07_R_data.txt", sep="\t", header=TRUE, comment.char="#",
                       na.strings=c("","NA"), stringsAsFactors=FALSE,
                       quote="", fill=FALSE)

Tte_data <- read.table("4_Tte_b3v07_R_data.txt", sep="\t", header=TRUE, comment.char="#",
                       na.strings=c("","NA"), stringsAsFactors=FALSE,
                       quote="", fill=FALSE)

Tpa_data <- read.table("5_Tpa_b3v07_R_data.txt", sep="\t", header=TRUE, comment.char="#",
                       na.strings=c("","NA"), stringsAsFactors=FALSE,
                       quote="", fill=FALSE)

Tge_data <- read.table("5_Tge_b3v07_R_data.txt", sep="\t", header=TRUE, comment.char="#",
                       na.strings=c("","NA"), stringsAsFactors=FALSE,
                       quote="", fill=FALSE)

#head(data)

## Change level to make it like RepeatMasker output

Tps_data$Type <- factor(Tps_data$Type, levels = c('SINE', 'LINE/R2', 'LINE/Jockey-I', 'LINE/R1', 'LINE/LOA', 'LINE/L2', 'LINE/CR1', 'LINE/RTE', 'LINE', 'LINE/L1', 'LTR', 'LTR/Gypsy', 'LTR/Copia', 'RC/Helitron', 'DNA/Transib', 'DNA/TcMar', 'DNA/Sola', 'DNA/PiggyBac', 'DNA/P', 'DNA', 'DNA/Maverick', 'DNA/Kolobok', 'DNA/hAT', 'DNA/Harbinger', 'DNA/Ginger', 'DNA/Academ', 'Other'))
Tdi_data$Type <- factor(Tdi_data$Type, levels = c('SINE', 'LINE/R2', 'LINE/Jockey-I', 'LINE/R1', 'LINE/LOA', 'LINE/L2', 'LINE/CR1', 'LINE/RTE', 'LINE', 'LINE/L1', 'LTR', 'LTR/Gypsy', 'LTR/Copia', 'RC/Helitron', 'DNA/Transib', 'DNA/TcMar', 'DNA/Sola', 'DNA/PiggyBac', 'DNA/P', 'DNA', 'DNA/Maverick', 'DNA/Kolobok', 'DNA/hAT', 'DNA/Harbinger', 'DNA/Ginger', 'DNA/Academ', 'Other'))
Tcm_data$Type <- factor(Tcm_data$Type, levels = c('SINE', 'LINE/R2', 'LINE/Jockey-I', 'LINE/R1', 'LINE/LOA', 'LINE/L2', 'LINE/CR1', 'LINE/RTE', 'LINE', 'LINE/L1', 'LTR', 'LTR/Gypsy', 'LTR/Copia', 'RC/Helitron', 'DNA/Transib', 'DNA/TcMar', 'DNA/Sola', 'DNA/PiggyBac', 'DNA/P', 'DNA', 'DNA/Maverick', 'DNA/Kolobok', 'DNA/hAT', 'DNA/Harbinger', 'DNA/Ginger', 'DNA/Academ', 'Other'))
Tsi_data$Type <- factor(Tsi_data$Type, levels = c('SINE', 'LINE/R2', 'LINE/Jockey-I', 'LINE/R1', 'LINE/LOA', 'LINE/L2', 'LINE/CR1', 'LINE/RTE', 'LINE', 'LINE/L1', 'LTR', 'LTR/Gypsy', 'LTR/Copia', 'RC/Helitron', 'DNA/Transib', 'DNA/TcMar', 'DNA/Sola', 'DNA/PiggyBac', 'DNA/P', 'DNA', 'DNA/Maverick', 'DNA/Kolobok', 'DNA/hAT', 'DNA/Harbinger', 'DNA/Ginger', 'DNA/Academ', 'Other'))
Tce_data$Type <- factor(Tce_data$Type, levels = c('SINE', 'LINE/R2', 'LINE/Jockey-I', 'LINE/R1', 'LINE/LOA', 'LINE/L2', 'LINE/CR1', 'LINE/RTE', 'LINE', 'LINE/L1', 'LTR', 'LTR/Gypsy', 'LTR/Copia', 'RC/Helitron', 'DNA/Transib', 'DNA/TcMar', 'DNA/Sola', 'DNA/PiggyBac', 'DNA/P', 'DNA', 'DNA/Maverick', 'DNA/Kolobok', 'DNA/hAT', 'DNA/Harbinger', 'DNA/Ginger', 'DNA/Academ', 'Other'))
Tms_data$Type <- factor(Tms_data$Type, levels = c('SINE', 'LINE/R2', 'LINE/Jockey-I', 'LINE/R1', 'LINE/LOA', 'LINE/L2', 'LINE/CR1', 'LINE/RTE', 'LINE', 'LINE/L1', 'LTR', 'LTR/Gypsy', 'LTR/Copia', 'RC/Helitron', 'DNA/Transib', 'DNA/TcMar', 'DNA/Sola', 'DNA/PiggyBac', 'DNA/P', 'DNA', 'DNA/Maverick', 'DNA/Kolobok', 'DNA/hAT', 'DNA/Harbinger', 'DNA/Ginger', 'DNA/Academ', 'Other'))
Tbi_data$Type <- factor(Tbi_data$Type, levels = c('SINE', 'LINE/R2', 'LINE/Jockey-I', 'LINE/R1', 'LINE/LOA', 'LINE/L2', 'LINE/CR1', 'LINE/RTE', 'LINE', 'LINE/L1', 'LTR', 'LTR/Gypsy', 'LTR/Copia', 'RC/Helitron', 'DNA/Transib', 'DNA/TcMar', 'DNA/Sola', 'DNA/PiggyBac', 'DNA/P', 'DNA', 'DNA/Maverick', 'DNA/Kolobok', 'DNA/hAT', 'DNA/Harbinger', 'DNA/Ginger', 'DNA/Academ', 'Other'))
Tte_data$Type <- factor(Tte_data$Type, levels = c('SINE', 'LINE/R2', 'LINE/Jockey-I', 'LINE/R1', 'LINE/LOA', 'LINE/L2', 'LINE/CR1', 'LINE/RTE', 'LINE', 'LINE/L1', 'LTR', 'LTR/Gypsy', 'LTR/Copia', 'RC/Helitron', 'DNA/Transib', 'DNA/TcMar', 'DNA/Sola', 'DNA/PiggyBac', 'DNA/P', 'DNA', 'DNA/Maverick', 'DNA/Kolobok', 'DNA/hAT', 'DNA/Harbinger', 'DNA/Ginger', 'DNA/Academ', 'Other'))
Tpa_data$Type <- factor(Tpa_data$Type, levels = c('SINE', 'LINE/R2', 'LINE/Jockey-I', 'LINE/R1', 'LINE/LOA', 'LINE/L2', 'LINE/CR1', 'LINE/RTE', 'LINE', 'LINE/L1', 'LTR', 'LTR/Gypsy', 'LTR/Copia', 'RC/Helitron', 'DNA/Transib', 'DNA/TcMar', 'DNA/Sola', 'DNA/PiggyBac', 'DNA/P', 'DNA', 'DNA/Maverick', 'DNA/Kolobok', 'DNA/hAT', 'DNA/Harbinger', 'DNA/Ginger', 'DNA/Academ', 'Other'))
Tge_data$Type <- factor(Tge_data$Type, levels = c('SINE', 'LINE/R2', 'LINE/Jockey-I', 'LINE/R1', 'LINE/LOA', 'LINE/L2', 'LINE/CR1', 'LINE/RTE', 'LINE', 'LINE/L1', 'LTR', 'LTR/Gypsy', 'LTR/Copia', 'RC/Helitron', 'DNA/Transib', 'DNA/TcMar', 'DNA/Sola', 'DNA/PiggyBac', 'DNA/P', 'DNA', 'DNA/Maverick', 'DNA/Kolobok', 'DNA/hAT', 'DNA/Harbinger', 'DNA/Ginger', 'DNA/Academ', 'Other'))

## Change color to make it like RepeatMasker output

fill <- c('#9F1FF0', '#A3C5DE', '#8FA1CF', '#848FC8', '#797EC0', '#625CB1', '#483AA2', '#38299A', '#251792', '#00008B', '#65BD61', '#489E42', '#3A8F33', '#FF00FF', '#FF9972', '#FF936C', '#FF8D65', '#FF865E', '#FF7F57', '#FF6A42', '#FF623B', '#FF5A34', '#FF512D', '#FF4825', '#FF3D1E', '#FF0000', '#4D4D4D')

# Stacked barplot

png(filename="Tps_data.png",res = 180, height = 977, width = 1877)
Tps_plot <- ggplot(Tps_data, aes(fill=Type, y=Genome_percent, x=Kimura)) + 
  geom_bar( stat="identity") + scale_fill_manual(values=fill) + 
  labs(x="Kimura substitution level (CpG adjusted)", y="Percent of genome") +
  coord_cartesian(ylim = c(0, 1.3)) + theme(legend.position = "none") +
  ggtitle("Tps")
Tps_plot
dev.off()

png(filename="Tdi_plot.png",res = 180, height = 977, width = 1877)
Tdi_plot <- ggplot(Tdi_data, aes(fill=Type, y=Genome_percent, x=Kimura)) + 
  geom_bar( stat="identity") + scale_fill_manual(values=fill) + 
  labs(x="Kimura substitution level (CpG adjusted)", y="Percent of genome") +
  coord_cartesian(ylim = c(0, 1.3)) + theme(legend.position = "none") +
  ggtitle("Tdi")
Tdi_plot
dev.off()

png(filename="Tcm_plot.png",res = 180, height = 977, width = 1877)
Tcm_plot <- ggplot(Tcm_data, aes(fill=Type, y=Genome_percent, x=Kimura)) + 
  geom_bar( stat="identity") + scale_fill_manual(values=fill) + 
  labs(x="Kimura substitution level (CpG adjusted)", y="Percent of genome") +
  coord_cartesian(ylim = c(0, 1.3)) + theme(legend.position = "none") +
  ggtitle("Tcm")
Tcm_plot
dev.off()

png(filename="Tsi_plot.png",res = 180, height = 977, width = 1877)
Tsi_plot <- ggplot(Tsi_data, aes(fill=Type, y=Genome_percent, x=Kimura)) + 
  geom_bar( stat="identity") + scale_fill_manual(values=fill) + 
  labs(x="Kimura substitution level (CpG adjusted)", y="Percent of genome") +
  coord_cartesian(ylim = c(0, 1.3)) + theme(legend.position = "none") +
  ggtitle("Tsi")
Tsi_plot
dev.off()

png(filename="Tce_plot.png",res = 180, height = 977, width = 1877)
Tce_plot <- ggplot(Tce_data, aes(fill=Type, y=Genome_percent, x=Kimura)) + 
  geom_bar( stat="identity") + scale_fill_manual(values=fill) + 
  labs(x="Kimura substitution level (CpG adjusted)", y="Percent of genome") +
  coord_cartesian(ylim = c(0, 1.3)) + theme(legend.position = "none") +
  ggtitle("Tce")
Tce_plot
dev.off()

png(filename="Tms_plot.png",res = 180, height = 977, width = 1877)
Tms_plot <- ggplot(Tms_data, aes(fill=Type, y=Genome_percent, x=Kimura)) + 
  geom_bar( stat="identity") + scale_fill_manual(values=fill) + 
  labs(x="Kimura substitution level (CpG adjusted)", y="Percent of genome") +
  coord_cartesian(ylim = c(0, 1.3)) + theme(legend.position = "none") +
  ggtitle("Tms")
Tms_plot
dev.off()

png(filename="Tbi_plot.png",res = 180, height = 977, width = 1877)
Tbi_plot <- ggplot(Tbi_data, aes(fill=Type, y=Genome_percent, x=Kimura)) + 
  geom_bar( stat="identity") + scale_fill_manual(values=fill) + 
  labs(x="Kimura substitution level (CpG adjusted)", y="Percent of genome") +
  coord_cartesian(ylim = c(0, 1.3)) + theme(legend.position = "none") +
  ggtitle("Tbi")
Tbi_plot
dev.off()

png(filename="Tte_plot.png",res = 180, height = 977, width = 1877)
Tte_plot <- ggplot(Tte_data, aes(fill=Type, y=Genome_percent, x=Kimura)) + 
  geom_bar( stat="identity") + scale_fill_manual(values=fill) + 
  labs(x="Kimura substitution level (CpG adjusted)", y="Percent of genome") +
  coord_cartesian(ylim = c(0, 1.3)) + theme(legend.position = "none") +
  ggtitle("Tte")
Tte_plot
dev.off()

png(filename="Tpa_plot.png",res = 180, height = 977, width = 1877)
Tpa_plot <- ggplot(Tpa_data, aes(fill=Type, y=Genome_percent, x=Kimura)) + 
  geom_bar( stat="identity") + scale_fill_manual(values=fill) + 
  labs(x="Kimura substitution level (CpG adjusted)", y="Percent of genome") +
  coord_cartesian(ylim = c(0, 1.3)) + theme(legend.position = "none") +
  ggtitle("Tpa")
Tpa_plot
dev.off()

png(filename="Tge_plot.png",res = 180, height = 977, width = 1877)
Tge_plot <- ggplot(Tge_data, aes(fill=Type, y=Genome_percent, x=Kimura)) + 
  geom_bar( stat="identity") + scale_fill_manual(values=fill) + 
  labs(x="Kimura substitution level (CpG adjusted)", y="Percent of genome") +
  coord_cartesian(ylim = c(0, 1.3)) + theme(legend.position = "none") +
  ggtitle("Tge")
Tge_plot
dev.off()

#ggarrange(Tps_plot, Tdi_plot, Tcm_plot, Tsi_plot, Tce_plot, Tms_plot, Tbi_plot, Tte_plot, Tpa_plot, Tge_plot, ncol=2, nrow=5, common.legend = TRUE, legend="bottom")

my_legend <- get_legend(Tpa_plot)
png(filename="legend.png",res = 180, height = 977, width = 1877)
as_ggplot(my_legend)
dev.off()

