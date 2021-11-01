library(ggplot2)
library(cowplot)
library(ggpubr)

#v1.b

# https://www.r-graph-gallery.com/48-grouped-barplot-with-ggplot2/

setwd("/Users/admin/Documents/schwander/timema/TE/landscapes_knownTEs")

Tps_data <- read.table("1_Tps_b3v08_R_data.txt", sep="\t", header=TRUE, comment.char="#",
                   na.strings=c("","NA"), stringsAsFactors=FALSE,
                   quote="", fill=FALSE)

Tdi_data <- read.table("1_Tdi_b3v08_R_data.txt", sep="\t", header=TRUE, comment.char="#",
                       na.strings=c("","NA"), stringsAsFactors=FALSE,
                       quote="", fill=FALSE)

Tcm_data <- read.table("2_Tcm_b3v08_R_data.txt", sep="\t", header=TRUE, comment.char="#",
                       na.strings=c("","NA"), stringsAsFactors=FALSE,
                       quote="", fill=FALSE)

Tsi_data <- read.table("2_Tsi_b3v08_R_data.txt", sep="\t", header=TRUE, comment.char="#",
                       na.strings=c("","NA"), stringsAsFactors=FALSE,
                       quote="", fill=FALSE)

Tce_data <- read.table("3_Tce_b3v08_R_data.txt", sep="\t", header=TRUE, comment.char="#",
                       na.strings=c("","NA"), stringsAsFactors=FALSE,
                       quote="", fill=FALSE)

Tms_data <- read.table("3_Tms_b3v08_R_data.txt", sep="\t", header=TRUE, comment.char="#",
                       na.strings=c("","NA"), stringsAsFactors=FALSE,
                       quote="", fill=FALSE)

Tbi_data <- read.table("4_Tbi_b3v08_R_data.txt", sep="\t", header=TRUE, comment.char="#",
                       na.strings=c("","NA"), stringsAsFactors=FALSE,
                       quote="", fill=FALSE)

Tte_data <- read.table("4_Tte_b3v08_R_data.txt", sep="\t", header=TRUE, comment.char="#",
                       na.strings=c("","NA"), stringsAsFactors=FALSE,
                       quote="", fill=FALSE)

Tpa_data <- read.table("5_Tpa_b3v08_R_data.txt", sep="\t", header=TRUE, comment.char="#",
                       na.strings=c("","NA"), stringsAsFactors=FALSE,
                       quote="", fill=FALSE)

Tge_data <- read.table("5_Tge_b3v08_R_data.txt", sep="\t", header=TRUE, comment.char="#",
                       na.strings=c("","NA"), stringsAsFactors=FALSE,
                       quote="", fill=FALSE)

#head(data)

## Change level to make it like RepeatMasker output

Tps_data$Type <- factor(Tps_data$Type, levels = c('SINE', 'LINE/Crack', 'LINE/Tx1', 'LINE/Poseidon', 'LINE/R2', 'LINE/Jockey-I', 'LINE/R1', 'LINE/LOA', 'LINE/L2', 'LINE/CR1', 'LINE/RTE', 'LINE', 'LINE/L1', 'LTR/Penelope', 'LTR', 'LTR/Gypsy', 'LTR/Copia', 'LTR/BEL', 'RC/Helitron', 'DNA/ISL2EU', 'DNA/MITE', 'DNA/MuDr', 'DNA/Transib', 'DNA/TcMar', 'DNA/Sola', 'DNA/PiggyBac', 'DNA/P', 'DNA', 'DNA/Maverick', 'DNA/Kolobok', 'DNA/hAT', 'DNA/Harbinger', 'DNA/Ginger', 'DNA/Chapaev', 'DNA/Academ', 'Other'))

Tdi_data$Type <- factor(Tdi_data$Type, levels = c('SINE', 'LINE/Crack', 'LINE/Tx1', 'LINE/Poseidon', 'LINE/R2', 'LINE/Jockey-I', 'LINE/R1', 'LINE/LOA', 'LINE/L2', 'LINE/CR1', 'LINE/RTE', 'LINE', 'LINE/L1', 'LTR/Penelope', 'LTR', 'LTR/Gypsy', 'LTR/Copia', 'LTR/BEL', 'RC/Helitron', 'DNA/ISL2EU', 'DNA/MITE', 'DNA/MuDr', 'DNA/Transib', 'DNA/TcMar', 'DNA/Sola', 'DNA/PiggyBac', 'DNA/P', 'DNA', 'DNA/Maverick', 'DNA/Kolobok', 'DNA/hAT', 'DNA/Harbinger', 'DNA/Ginger', 'DNA/Chapaev', 'DNA/Academ', 'Other'))

Tcm_data$Type <- factor(Tcm_data$Type, levels = c('SINE', 'LINE/Crack', 'LINE/Tx1', 'LINE/Poseidon', 'LINE/R2', 'LINE/Jockey-I', 'LINE/R1', 'LINE/LOA', 'LINE/L2', 'LINE/CR1', 'LINE/RTE', 'LINE', 'LINE/L1', 'LTR/Penelope', 'LTR', 'LTR/Gypsy', 'LTR/Copia', 'LTR/BEL', 'RC/Helitron', 'DNA/ISL2EU', 'DNA/MITE', 'DNA/MuDr', 'DNA/Transib', 'DNA/TcMar', 'DNA/Sola', 'DNA/PiggyBac', 'DNA/P', 'DNA', 'DNA/Maverick', 'DNA/Kolobok', 'DNA/hAT', 'DNA/Harbinger', 'DNA/Ginger', 'DNA/Chapaev', 'DNA/Academ', 'Other'))

Tsi_data$Type <- factor(Tsi_data$Type, levels = c('SINE', 'LINE/Crack', 'LINE/Tx1', 'LINE/Poseidon', 'LINE/R2', 'LINE/Jockey-I', 'LINE/R1', 'LINE/LOA', 'LINE/L2', 'LINE/CR1', 'LINE/RTE', 'LINE', 'LINE/L1', 'LTR/Penelope', 'LTR', 'LTR/Gypsy', 'LTR/Copia', 'LTR/BEL', 'RC/Helitron', 'DNA/ISL2EU', 'DNA/MITE', 'DNA/MuDr', 'DNA/Transib', 'DNA/TcMar', 'DNA/Sola', 'DNA/PiggyBac', 'DNA/P', 'DNA', 'DNA/Maverick', 'DNA/Kolobok', 'DNA/hAT', 'DNA/Harbinger', 'DNA/Ginger', 'DNA/Chapaev', 'DNA/Academ', 'Other'))

Tce_data$Type <- factor(Tce_data$Type, levels = c('SINE', 'LINE/Crack', 'LINE/Tx1', 'LINE/Poseidon', 'LINE/R2', 'LINE/Jockey-I', 'LINE/R1', 'LINE/LOA', 'LINE/L2', 'LINE/CR1', 'LINE/RTE', 'LINE', 'LINE/L1', 'LTR/Penelope', 'LTR', 'LTR/Gypsy', 'LTR/Copia', 'LTR/BEL', 'RC/Helitron', 'DNA/ISL2EU', 'DNA/MITE', 'DNA/MuDr', 'DNA/Transib', 'DNA/TcMar', 'DNA/Sola', 'DNA/PiggyBac', 'DNA/P', 'DNA', 'DNA/Maverick', 'DNA/Kolobok', 'DNA/hAT', 'DNA/Harbinger', 'DNA/Ginger', 'DNA/Chapaev', 'DNA/Academ', 'Other'))

Tms_data$Type <- factor(Tms_data$Type, levels = c('SINE', 'LINE/Crack', 'LINE/Tx1', 'LINE/Poseidon', 'LINE/R2', 'LINE/Jockey-I', 'LINE/R1', 'LINE/LOA', 'LINE/L2', 'LINE/CR1', 'LINE/RTE', 'LINE', 'LINE/L1', 'LTR/Penelope', 'LTR', 'LTR/Gypsy', 'LTR/Copia', 'LTR/BEL', 'RC/Helitron', 'DNA/ISL2EU', 'DNA/MITE', 'DNA/MuDr', 'DNA/Transib', 'DNA/TcMar', 'DNA/Sola', 'DNA/PiggyBac', 'DNA/P', 'DNA', 'DNA/Maverick', 'DNA/Kolobok', 'DNA/hAT', 'DNA/Harbinger', 'DNA/Ginger', 'DNA/Chapaev', 'DNA/Academ', 'Other'))

Tbi_data$Type <- factor(Tbi_data$Type, levels = c('SINE', 'LINE/Crack', 'LINE/Tx1', 'LINE/Poseidon', 'LINE/R2', 'LINE/Jockey-I', 'LINE/R1', 'LINE/LOA', 'LINE/L2', 'LINE/CR1', 'LINE/RTE', 'LINE', 'LINE/L1', 'LTR/Penelope', 'LTR', 'LTR/Gypsy', 'LTR/Copia', 'LTR/BEL', 'RC/Helitron', 'DNA/ISL2EU', 'DNA/MITE', 'DNA/MuDr', 'DNA/Transib', 'DNA/TcMar', 'DNA/Sola', 'DNA/PiggyBac', 'DNA/P', 'DNA', 'DNA/Maverick', 'DNA/Kolobok', 'DNA/hAT', 'DNA/Harbinger', 'DNA/Ginger', 'DNA/Chapaev', 'DNA/Academ', 'Other'))

Tte_data$Type <- factor(Tte_data$Type, levels = c('SINE', 'LINE/Crack', 'LINE/Tx1', 'LINE/Poseidon', 'LINE/R2', 'LINE/Jockey-I', 'LINE/R1', 'LINE/LOA', 'LINE/L2', 'LINE/CR1', 'LINE/RTE', 'LINE', 'LINE/L1', 'LTR/Penelope', 'LTR', 'LTR/Gypsy', 'LTR/Copia', 'LTR/BEL', 'RC/Helitron', 'DNA/ISL2EU', 'DNA/MITE', 'DNA/MuDr', 'DNA/Transib', 'DNA/TcMar', 'DNA/Sola', 'DNA/PiggyBac', 'DNA/P', 'DNA', 'DNA/Maverick', 'DNA/Kolobok', 'DNA/hAT', 'DNA/Harbinger', 'DNA/Ginger', 'DNA/Chapaev', 'DNA/Academ', 'Other'))

Tpa_data$Type <- factor(Tpa_data$Type, levels = c('SINE', 'LINE/Crack', 'LINE/Tx1', 'LINE/Poseidon', 'LINE/R2', 'LINE/Jockey-I', 'LINE/R1', 'LINE/LOA', 'LINE/L2', 'LINE/CR1', 'LINE/RTE', 'LINE', 'LINE/L1', 'LTR/Penelope', 'LTR', 'LTR/Gypsy', 'LTR/Copia', 'LTR/BEL', 'RC/Helitron', 'DNA/ISL2EU', 'DNA/MITE', 'DNA/MuDr', 'DNA/Transib', 'DNA/TcMar', 'DNA/Sola', 'DNA/PiggyBac', 'DNA/P', 'DNA', 'DNA/Maverick', 'DNA/Kolobok', 'DNA/hAT', 'DNA/Harbinger', 'DNA/Ginger', 'DNA/Chapaev', 'DNA/Academ', 'Other'))

Tge_data$Type <- factor(Tge_data$Type, levels = c('SINE', 'LINE/Crack', 'LINE/Tx1', 'LINE/Poseidon', 'LINE/R2', 'LINE/Jockey-I', 'LINE/R1', 'LINE/LOA', 'LINE/L2', 'LINE/CR1', 'LINE/RTE', 'LINE', 'LINE/L1', 'LTR/Penelope', 'LTR', 'LTR/Gypsy', 'LTR/Copia', 'LTR/BEL', 'RC/Helitron', 'DNA/ISL2EU', 'DNA/MITE', 'DNA/MuDr', 'DNA/Transib', 'DNA/TcMar', 'DNA/Sola', 'DNA/PiggyBac', 'DNA/P', 'DNA', 'DNA/Maverick', 'DNA/Kolobok', 'DNA/hAT', 'DNA/Harbinger', 'DNA/Ginger', 'DNA/Chapaev', 'DNA/Academ', 'Other'))

## Change color to make it like RepeatMasker output


fill <- c('#9F1FF0', '#c5d4eb', '#B5C9E8', '#C1DFE8', '#A3C5DE', '#8FA1CF', '#848FC8', '#797EC0', '#625CB1', '#483AA2', '#38299A', '#251792', '#00008B', '#57A157', '#65BD61', '#489E42', '#3A8F33', '#013601', '#FF00FF', '#e3baaa', '#D19A84', '#fad6c8', '#FF9972', '#FF936C', '#FF8D65', '#FF865E', '#FF7F57', '#FF6A42', '#FF623B', '#FF5A34', '#FF512D', '#FF4825', '#FF3D1E', '#FF2B0F', '#FF0000', '#4D4D4D')

# Stacked barplot

png(filename="Tps_data.png",res = 180, height = 977, width = 1877)
Tps_plot <- ggplot(Tps_data, aes(fill=Type, y=Genome_percent, x=Kimura)) + 
  geom_bar( stat="identity") + scale_fill_manual(values=fill) + 
  labs(x="Kimura substitution level (CpG adjusted)", y="Percent of genome") +
  coord_cartesian(ylim = c(0, 2)) +
  scale_y_continuous(breaks = seq(0, 2, 0.2)) +
  labs(title = "T. poppensis") +
  theme_bw() + theme(legend.position = "none", plot.title = element_text(color = "red", size = 15, face = "italic"))
Tps_plot
dev.off()

png(filename="Tdi_plot.png",res = 180, height = 977, width = 1877)
Tdi_plot <- ggplot(Tdi_data, aes(fill=Type, y=Genome_percent, x=Kimura)) + 
  geom_bar( stat="identity") + scale_fill_manual(values=fill) + 
  labs(x="Kimura substitution level (CpG adjusted)", y="Percent of genome") +
  coord_cartesian(ylim = c(0, 2)) +
  scale_y_continuous(breaks = seq(0, 2, 0.2)) +
  labs(title = "T. douglasi") +
  theme_bw() + theme(legend.position = "none", plot.title = element_text(color = "blue", size = 15, face = "italic"))
Tdi_plot
dev.off()

png(filename="Tcm_plot.png",res = 180, height = 977, width = 1877)
Tcm_plot <- ggplot(Tcm_data, aes(fill=Type, y=Genome_percent, x=Kimura)) + 
  geom_bar( stat="identity") + scale_fill_manual(values=fill) + 
  labs(x="Kimura substitution level (CpG adjusted)", y="Percent of genome") +
  coord_cartesian(ylim = c(0, 2)) +
  scale_y_continuous(breaks = seq(0, 2, 0.2)) +
  labs(title = "T. californicum") +
  theme_bw() + theme(legend.position = "none", plot.title = element_text(color = "red", size = 15, face = "italic"))
Tcm_plot
dev.off()

png(filename="Tsi_plot.png",res = 180, height = 977, width = 1877)
Tsi_plot <- ggplot(Tsi_data, aes(fill=Type, y=Genome_percent, x=Kimura)) + 
  geom_bar( stat="identity") + scale_fill_manual(values=fill) + 
  labs(x="Kimura substitution level (CpG adjusted)", y="Percent of genome") +
  coord_cartesian(ylim = c(0, 2)) +
  scale_y_continuous(breaks = seq(0, 2, 0.2)) +
  labs(title = "T. shepardi") +
  theme_bw() + theme(legend.position = "none", plot.title = element_text(color = "blue", size = 15, face = "italic"))
Tsi_plot
dev.off()

png(filename="Tce_plot.png",res = 180, height = 977, width = 1877)
Tce_plot <- ggplot(Tce_data, aes(fill=Type, y=Genome_percent, x=Kimura)) + 
  geom_bar( stat="identity") + scale_fill_manual(values=fill) + 
  labs(x="Kimura substitution level (CpG adjusted)", y="Percent of genome") +
  coord_cartesian(ylim = c(0, 2)) +
  scale_y_continuous(breaks = seq(0, 2, 0.2)) +
  labs(title = "T. cristinae") +
  theme_bw() + theme(legend.position = "none", plot.title = element_text(color = "red", size = 15, face = "italic"))
Tce_plot
dev.off()

png(filename="Tms_plot.png",res = 180, height = 977, width = 1877)
Tms_plot <- ggplot(Tms_data, aes(fill=Type, y=Genome_percent, x=Kimura)) + 
  geom_bar( stat="identity") + scale_fill_manual(values=fill) + 
  labs(x="Kimura substitution level (CpG adjusted)", y="Percent of genome") +
  coord_cartesian(ylim = c(0, 2)) +
  scale_y_continuous(breaks = seq(0, 2, 0.2)) +
  labs(title = "T. monikensis") +
  theme_bw() + theme(legend.position = "none", plot.title = element_text(color = "blue", size = 15, face = "italic"))
Tms_plot
dev.off()

png(filename="Tbi_plot.png",res = 180, height = 977, width = 1877)
Tbi_plot <- ggplot(Tbi_data, aes(fill=Type, y=Genome_percent, x=Kimura)) + 
  geom_bar( stat="identity") + scale_fill_manual(values=fill) + 
  labs(x="Kimura substitution level (CpG adjusted)", y="Percent of genome") +
  coord_cartesian(ylim = c(0, 2)) +
  scale_y_continuous(breaks = seq(0, 2, 0.2)) +
  labs(title = "T. bartmani") +
  theme_bw() + theme(legend.position = "none", plot.title = element_text(color = "red", size = 15, face = "italic"))
Tbi_plot
dev.off()

png(filename="Tte_plot.png",res = 180, height = 977, width = 1877)
Tte_plot <- ggplot(Tte_data, aes(fill=Type, y=Genome_percent, x=Kimura)) + 
  geom_bar( stat="identity") + scale_fill_manual(values=fill) + 
  labs(x="Kimura substitution level (CpG adjusted)", y="Percent of genome") +
  coord_cartesian(ylim = c(0, 2)) +
  scale_y_continuous(breaks = seq(0, 2, 0.2)) +
  labs(title = "T. tahoe") +
  theme_bw() + theme(legend.position = "none", plot.title = element_text(color = "blue", size = 15, face = "italic"))
Tte_plot
dev.off()

png(filename="Tpa_plot.png",res = 180, height = 977, width = 1877)
Tpa_plot <- ggplot(Tpa_data, aes(fill=Type, y=Genome_percent, x=Kimura)) + 
  geom_bar( stat="identity") + scale_fill_manual(values=fill) + 
  labs(x="Kimura substitution level (CpG adjusted)", y="Percent of genome") +
  coord_cartesian(ylim = c(0, 2)) +
  scale_y_continuous(breaks = seq(0, 2, 0.2)) +
  labs(title = "T. podura") +
  theme_bw() + theme(legend.position = "none", plot.title = element_text(color = "red", size = 15, face = "italic"))
Tpa_plot
dev.off()

png(filename="Tge_plot.png",res = 180, height = 977, width = 1877)
Tge_plot <- ggplot(Tge_data, aes(fill=Type, y=Genome_percent, x=Kimura)) + 
  geom_bar( stat="identity") + scale_fill_manual(values=fill) + 
  labs(x="Kimura substitution level (CpG adjusted)", y="Percent of genome") +
  coord_cartesian(ylim = c(0, 2)) +
  scale_y_continuous(breaks = seq(0, 2, 0.2)) +
  labs(title = "T. genevievae") +
  theme_bw() + theme(legend.position = "none", plot.title = element_text(color = "blue", size = 15, face = "italic"))
Tge_plot
dev.off()


p1 <- ggarrange(Tps_plot, Tdi_plot, ncol=2) 
p2 <- ggarrange(Tcm_plot, Tsi_plot, ncol=2) 
p3 <- ggarrange(Tce_plot, Tms_plot, ncol=2) 
p4 <- ggarrange(Tbi_plot, Tte_plot, ncol=2) 
p5 <- ggarrange(Tpa_plot, Tge_plot, ncol=2, common.legend = TRUE, legend = "bottom") 


pdf("TE_landscapes_all.pdf", width = 9, height = 15)
plot_grid(p1, p2, p3, p4, p5, nrow = 5,  rel_heights = c(1,1,1,1,1.8))
dev.off()

png(filename = "TE_landscapes_all.png", width = 9, height = 15, units = "in", bg = "white", res = 600)
plot_grid(p1, p2, p3, p4, p5, nrow = 5,  rel_heights = c(1,1,1,1,1.8))
dev.off()


pdf("TE_landscapes_all_2.pdf", width = 9, height = 15)
plot_grid(p1, p2, p3, p4, p5, nrow = 5,  rel_heights = c(1,1,1,1,1.8))
dev.off()

png(filename = "TE_landscapes_all_2.png", width = 8, height = 12, units = "in", bg = "white", res = 600)
plot_grid(p1, p2, p3, p4, p5, nrow = 5,  rel_heights = c(1,1,1,1,2.05))
dev.off()