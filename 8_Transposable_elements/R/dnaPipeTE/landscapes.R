library(ggplot2)
library(ggpubr)


setwd("/Users/admin/Documents/schwander/timema/TE/dnaPipeTE")

# Read the corresponding factor order and color table.

colors_file = "1_Tdi_factors_and_colors.txt"
fac_col = read.table(colors_file)

# Except for Tsi that miss LINE/L1

Tsi_colors_file = "2_Tsi_factors_and_colors.txt"
Tsi_fac_col = read.table(Tsi_colors_file)

#  head(fac_col)
#V1      V2      V3
#1    DNA/Academ ClassII #FF0000
#2           DNA ClassII #FF6A42

# Read dnaPipeTE landscape results

Tps_file = "1_Tps_reads_landscape.txt"
Tdi_file = "1_Tdi_reads_landscape.txt"
Tcm_file = "2_Tcm_reads_landscape.txt"
Tsi_file = "2_Tsi_reads_landscape.txt"
Tce_file = "3_Tce_reads_landscape.txt"
Tms_file = "3_Tms_reads_landscape.txt"
Tbi_file = "4_Tbi_reads_landscape.txt"
Tte_file = "4_Tte_reads_landscape.txt"
Tpa_file = "5_Tpa_reads_landscape.txt"
Tge_file = "5_Tge_reads_landscape.txt"

Tps_land = read.table(Tps_file)
Tdi_land = read.table(Tdi_file)
Tcm_land = read.table(Tcm_file)
Tsi_land = read.table(Tsi_file)
Tce_land = read.table(Tce_file)
Tms_land = read.table(Tms_file)
Tbi_land = read.table(Tbi_file)
Tte_land = read.table(Tte_file)
Tpa_land = read.table(Tpa_file)
Tge_land = read.table(Tge_file)

# Rename columns and create divergence information

names(Tps_land)=c("id", "annot", "fam1", "fam")
names(Tdi_land)=c("id", "annot", "fam1", "fam")
names(Tcm_land)=c("id", "annot", "fam1", "fam")
names(Tsi_land)=c("id", "annot", "fam1", "fam")
names(Tce_land)=c("id", "annot", "fam1", "fam")
names(Tms_land)=c("id", "annot", "fam1", "fam")
names(Tbi_land)=c("id", "annot", "fam1", "fam")
names(Tte_land)=c("id", "annot", "fam1", "fam")
names(Tpa_land)=c("id", "annot", "fam1", "fam")
names(Tge_land)=c("id", "annot", "fam1", "fam")

Tps_land$div = 100 - Tps_land$id
Tdi_land$div = 100 - Tdi_land$id
Tcm_land$div = 100 - Tcm_land$id
Tsi_land$div = 100 - Tsi_land$id
Tce_land$div = 100 - Tce_land$id
Tms_land$div = 100 - Tms_land$id
Tbi_land$div = 100 - Tbi_land$id
Tte_land$div = 100 - Tte_land$id
Tpa_land$div = 100 - Tpa_land$id
Tge_land$div = 100 - Tge_land$id

# Select only type of TEs that are in colors_file, the others are transformed
#to <NA>.
Tps_land$fam1 = factor(Tps_land$fam1, levels = as.character(fac_col$V1))
Tdi_land$fam1 = factor(Tdi_land$fam1, levels = as.character(fac_col$V1))
Tcm_land$fam1 = factor(Tcm_land$fam1, levels = as.character(fac_col$V1))
Tsi_land$fam1 = factor(Tsi_land$fam1, levels = as.character(Tsi_fac_col$V1))
Tce_land$fam1 = factor(Tce_land$fam1, levels = as.character(fac_col$V1))
Tms_land$fam1 = factor(Tms_land$fam1, levels = as.character(fac_col$V1))
Tbi_land$fam1 = factor(Tbi_land$fam1, levels = as.character(fac_col$V1))
Tte_land$fam1 = factor(Tte_land$fam1, levels = as.character(fac_col$V1))
Tpa_land$fam1 = factor(Tpa_land$fam1, levels = as.character(fac_col$V1))
Tge_land$fam1 = factor(Tge_land$fam1, levels = as.character(fac_col$V1))

#> head(Tdi_land)
#id annot    fam1 fam div
#1 100   DTA DNA/hAT DNA   0
#2 100   DTA DNA/hAT DNA   0
#..
#178  93.75   DXX         <NA>          DNA  6.25

# Plot the landscape graph

png(filename="Tps_plot.png",res = 180, height = 1200, width = 600)
Tps_plot <- ggplot(Tps_land, aes(div, fill=fam1)) + geom_histogram(binwidth=1.1) + 
  scale_fill_manual(values=as.character(fac_col$V3) ,na.value="#616161") +
  labs(x="Divergence", y="Reads count") + scale_x_continuous(limits = c(0, 35)) +
  coord_cartesian(ylim = c(0, 200000)) + theme(legend.position = "none") +
  ggtitle("Tps") 
Tps_plot
dev.off()

png(filename="Tdi_plot.png",res = 180, height = 1200, width = 600)
Tdi_plot <- ggplot(Tdi_land, aes(div, fill=fam1)) + geom_histogram(binwidth=1.1) + 
  scale_fill_manual(values=as.character(fac_col$V3) ,na.value="#616161") +
  labs(x="Divergence", y="Reads count") + scale_x_continuous(limits = c(0, 35)) +
  coord_cartesian(ylim = c(0, 200000)) + theme(legend.position = "none") +
  ggtitle("Tdi") 
Tdi_plot
dev.off()

png(filename="Tcm_plot.png",res = 180, height = 1200, width = 600)
Tcm_plot <- ggplot(Tcm_land, aes(div, fill=fam1)) + geom_histogram(binwidth=1.1) + 
  scale_fill_manual(values=as.character(fac_col$V3) ,na.value="#616161") +
  labs(x="Divergence", y="Reads count") + scale_x_continuous(limits = c(0, 35)) +
  coord_cartesian(ylim = c(0, 200000)) + theme(legend.position = "none") +
  ggtitle("Tcm") 
Tcm_plot
dev.off()

png(filename="Tsi_plot.png",res = 180, height = 1200, width = 600)
Tsi_plot <- ggplot(Tsi_land, aes(div, fill=fam1)) + geom_histogram(binwidth=1.1) + 
  scale_fill_manual(values=as.character(Tsi_fac_col$V3) ,na.value="#616161") +
  labs(x="Divergence", y="Reads count") + scale_x_continuous(limits = c(0, 35)) +
  coord_cartesian(ylim = c(0, 200000)) + theme(legend.position = "none") +
  ggtitle("Tsi") 
Tsi_plot
dev.off()

png(filename="Tce_plot.png",res = 180, height = 1200, width = 600)
Tce_plot <- ggplot(Tce_land, aes(div, fill=fam1)) + geom_histogram(binwidth=1.1) + 
  scale_fill_manual(values=as.character(fac_col$V3) ,na.value="#616161") +
  labs(x="Divergence", y="Reads count") + scale_x_continuous(limits = c(0, 35)) +
  coord_cartesian(ylim = c(0, 200000)) + theme(legend.position = "none") +
  ggtitle("Tce") 
Tce_plot
dev.off()

png(filename="Tms_plot.png",res = 180, height = 1200, width = 600)
Tms_plot <- ggplot(Tms_land, aes(div, fill=fam1)) + geom_histogram(binwidth=1.1) + 
  scale_fill_manual(values=as.character(fac_col$V3) ,na.value="#616161") +
  labs(x="Divergence", y="Reads count") + scale_x_continuous(limits = c(0, 35)) +
  coord_cartesian(ylim = c(0, 200000)) + theme(legend.position = "none") +
  ggtitle("Tms") 
Tms_plot
dev.off()

png(filename="Tbi_plot.png",res = 180, height = 1200, width = 600)
Tbi_plot <- ggplot(Tbi_land, aes(div, fill=fam1)) + geom_histogram(binwidth=1.1) + 
  scale_fill_manual(values=as.character(fac_col$V3) ,na.value="#616161") +
  labs(x="Divergence", y="Reads count") + scale_x_continuous(limits = c(0, 35)) +
  coord_cartesian(ylim = c(0, 200000)) + theme(legend.position = "none") +
  ggtitle("Tbi") 
Tbi_plot
dev.off()

png(filename="Tte_plot.png",res = 180, height = 1200, width = 600)
Tte_plot <- ggplot(Tte_land, aes(div, fill=fam1)) + geom_histogram(binwidth=1.1) + 
  scale_fill_manual(values=as.character(fac_col$V3) ,na.value="#616161") +
  labs(x="Divergence", y="Reads count") + scale_x_continuous(limits = c(0, 35)) +
  coord_cartesian(ylim = c(0, 200000)) + theme(legend.position = "none") +
  ggtitle("Tte") 
Tte_plot
dev.off()

png(filename="Tpa_plot.png",res = 180, height = 1200, width = 600)
Tpa_plot <- ggplot(Tpa_land, aes(div, fill=fam1)) + geom_histogram(binwidth=1.1) + 
  scale_fill_manual(values=as.character(fac_col$V3) ,na.value="#616161") +
  labs(x="Divergence", y="Reads count") + scale_x_continuous(limits = c(0, 35)) +
  coord_cartesian(ylim = c(0, 200000)) + theme(legend.position = "none") +
  ggtitle("Tpa") 
Tpa_plot
dev.off()

png(filename="Tge_plot.png",res = 180, height = 1200, width = 600)
Tge_plot <- ggplot(Tge_land, aes(div, fill=fam1)) + geom_histogram(binwidth=1.1) + 
  scale_fill_manual(values=as.character(fac_col$V3) ,na.value="#616161") +
  labs(x="Divergence", y="Reads count") + scale_x_continuous(limits = c(0, 35)) +
  coord_cartesian(ylim = c(0, 200000)) + theme(legend.position = "none") +
  ggtitle("Tge") 
Tge_plot
dev.off()

p1 <- ggarrange(Tps_plot, Tdi_plot, ncol=2) 
p2 <- ggarrange(Tcm_plot, Tsi_plot, ncol=2) 
p3 <- ggarrange(Tce_plot, Tms_plot, ncol=2) 
p4 <- ggarrange(Tbi_plot, Tte_plot, ncol=2) 
p5 <- ggarrange(Tpa_plot, Tge_plot, ncol=2, common.legend = TRUE, legend = "bottom") 

pdf("timema_dnaPipeTE_landscape.pdf", width=10, height=30 )
p <- ggarrange(p1, p2, p3, p4, p5, nrow=5, common.legend = TRUE, legend="bottom") 
p
dev.off()


# orginal plot
ggplot(Tdi_land, aes(div, fill=fam1))+geom_histogram(binwidth=1.1)+
  labs(list(x="div", y="reads count"))+scale_x_continuous(limits = c(0, 35))+
  scale_fill_manual(values=as.character(fac_col$V3))+
  guides(fill=guide_legend(ncol=3))+
  theme(legend.direction ="vertical",legend.position = "bottom")
ggsave("landscape.pdf", height=12, width=6)



