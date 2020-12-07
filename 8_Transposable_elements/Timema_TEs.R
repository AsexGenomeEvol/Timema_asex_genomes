#script for Sinergia TE mapped content

library(ggplot2)
library(plyr)
library(scales)
library(lattice)
library(reshape2)
library(gridExtra)
library(grid)
library(vegan)
library(Hmisc)



###########
#FUNCTIONS#
###########


#permutation ANOVA

anovsF3<-function(reponse,v1,v2,v3,nbrep=5000){
  rep<-reponse
  var1<-v1
  var2<-v2
  var3<-v3
  data<-data.frame(rep=rep,var1=var1,var2=var2,var3=var3)
  obs<-data.frame(t(anova(lm(rep~var1+var2*var3,data=data))$F))
  names(obs)<-c("Factor 1","Factor 2","Factor 3","Interaction","")
  for (i in 2:nbrep) {
    obs[i,]<-t(anova(lm(sample(rep)~var1+var2*var3,data=data))$F)
  }
  pval<-NULL
  for (j in 1:4) {
    pval[j]<-sum(obs[,j]>=obs[1,j])/length(obs[,j])
  }
  par(mfrow=c(4,1))
  hist(obs$"Factor 1",xlab="F");abline(v=obs[1,1],lwd=2)
  hist(obs$"Factor 2",xlab="F");abline(v=obs[1,2],lwd=2)
  hist(obs$"Factor 3",xlab="F");abline(v=obs[1,3],lwd=2)
  hist(obs$"Interaction",xlab="F");abline(v=obs[1,4],lwd=2)
  
  return(pval)
}


#median and CI bootstrap for plot

median_cl_boot <- function(x, conf = 0.95) {
  lconf <- (1 - conf)/2
  uconf <- 1 - lconf
  require(boot)
  bmedian <- function(x, ind) median(x[ind])
  bt <- boot(x, bmedian, 5000)
  bb <- boot.ci(bt, type = "perc")
  data.frame(y = median(x), ymin = quantile(bt$t, lconf), ymax = quantile(bt$t, uconf))
}


########################
#Get data and formating#
########################
#old poly (non-merged)
#TEfreq = read.table(file="~/Dropbox/UNIL/Timema/TEs/mappingRM/polymorphism_reseq/count/All_unq.count.comb",sep='\t', header=T)

#genome version 8
#with merged TElibs genome only
#TEfreq = read.table(file="~/Dropbox/UNIL/Timema/TEs/mappingRM/new_merged/All_unq.count",sep='\t', header=T)

#All data with genome and reseq reads with merged TE lib and version 08 genome
TEfreq = read.table(file="~/Dropbox/UNIL/Timema/TEs/mappingRM/new_merged/All_unq.count.comb2",sep='\t', header=T)
#map_stats = read.table(file="~/Dropbox/UNIL/Timema/TEs/mappingRM/new_merged/samplenames_map_stats_pol",sep='\t', header=F)


#specify names and mode of repr
names_asex <- c("Tdi","Tsi","Tms","Tge","Tte")
match_asex <- paste(names_asex, collapse="|")
names_sex <- c("Tps","Tcm","Tce","Tpa","Tbi")
match_sex <- paste(names_sex, collapse="|")

#add repr mode info in new column
TEfreq$mode <- ifelse(grepl(match_asex,TEfreq$species), "asex",
                           ifelse(grepl(match_sex,TEfreq$species), "sex","Other"))
#add pair info
TEfreq$pair <- ifelse(grepl("Tdi|Tps",TEfreq$species), "Tdi-Tps",
                           ifelse(grepl("Tsi|Tcm",TEfreq$species), "Tsi-Tcm",
                                  ifelse(grepl("Tms|Tce",TEfreq$species), "Tms-Tce",
                                         ifelse(grepl("Tge|Tpa",TEfreq$species), "Tge-Tpa",
                                                ifelse(grepl("Tte|Tbi",TEfreq$species), "Tte-Tbi", "Other")))))

#modify data type
TEfreq$pair <- as.factor(TEfreq$pair)
TEfreq$mode <- as.factor(TEfreq$mode)
#TEfreq$frac <- as.numeric(TEfreq$frac)
TEfreq$repl <- as.factor(TEfreq$repl)

#colors
palette = read.table(file="~/Dropbox/UNIL/Timema/stats_new/pieColorTEs", sep="\t")
names(palette)=c("class", "col")

#change order of TE type to adjust plot (according to color-palette)
#levels(TEfreq$type)
#TEfreq <- arrange(TEfreq, type)

#exclude some classes or the genome
#TEfreq <- subset(TEfreq, type!= "Genome")
#TEfreq <- subset(TEfreq, type!= "Unknown-Unknown")
TEfreq <- subset(TEfreq, TEfam!= "Unknown")
TEfreq <- subset(TEfreq, TEfam!= "XXX") #exclude SINE and Other
TEfreq <- subset(TEfreq, TEfam!= "SSR") #exclude simple repeats
TEfreq <- subset(TEfreq, TEfam!= "HG")  #exclude potential host genes
#TEfreq <- subset(TEfreq, TEfam!= "RSX") #exclude SINE


#exclude the "bad" Tsi samples
TEfreq <- subset(TEfreq, species!= "2_Tsi" | repl!="1") #exclude 2_Tsi repl 1
TEfreq <- subset(TEfreq, species!= "2_Tsi" | repl!="2") #exclude 2_Tsi repl 2



#bin unclear TE superfamily assignment
#CAUTION!!! Either replace here for barplot or below for the oder overview, but not both!!!
#TEfreq$TEfam <- sapply(TEfreq$TEfam, function(x){gsub(pattern = "DTX|DXX|RIX|RLX|RXX", replacement = "unclear", x)}) #unclear TE superfamily assignment
#TEfreq$TEfam <- sapply(TEfreq$TEfam, function(x){gsub(pattern = "DXX|RXX", replacement = "unclear", x)}) #unclear order

#TEfreq$TEfam <- as.factor(TEfreq$TEfam)



###########
#SUMMARIES#
###########

#sum count for the species and replicates levels (summing total TE load per replicate and species)
TEfreq_agg <- aggregate(count~species+repl, TEfreq, sum)
Mapped_agg <- aggregate(mapped~species+repl, TEfreq, mean)
TEfreq_agg$mapped <- as.integer(Mapped_agg$mapped)
TEfreq_agg$frac <- TEfreq_agg$count / TEfreq_agg$mapped
#str(TEfreq_agg)


#add repr mode info in new column
TEfreq_agg$mode <- ifelse(grepl(match_asex,TEfreq_agg$species), "asex",
                      ifelse(grepl(match_sex,TEfreq_agg$species), "sex","Other"))
#add pair info
TEfreq_agg$pair <- ifelse(grepl("Tdi|Tps",TEfreq_agg$species), "Tdi-Tps",
                      ifelse(grepl("Tsi|Tcm",TEfreq_agg$species), "Tsi-Tcm",
                             ifelse(grepl("Tms|Tce",TEfreq_agg$species), "Tms-Tce",
                                    ifelse(grepl("Tge|Tpa",TEfreq_agg$species), "Tge-Tpa",
                                           ifelse(grepl("Tte|Tbi",TEfreq_agg$species), "Tte-Tbi", "Other")))))

#add Clade info
TEfreq_agg$clade <- ifelse(grepl("Tdi|Tps|Tsi|Tcm|Tms|Tce",TEfreq_agg$species), "C1",
                          ifelse(grepl("Tge|Tpa|Tte|Tbi",TEfreq_agg$species), "C2","Other"))


#modify data type
TEfreq_agg$pair <- as.factor(TEfreq_agg$pair)
TEfreq_agg$mode <- as.factor(TEfreq_agg$mode)
#TEfreq_agg$frac <- as.numeric(TEfreq_agg$frac)
TEfreq_agg$repl <- as.factor(TEfreq_agg$repl)
TEfreq_agg$clade <- as.factor(TEfreq_agg$clade)

#difference in TE load between clades 1 and 2
tapply(TEfreq_agg$frac, TEfreq_agg$clade, mean)




#######
#PLOTS#
#######


#plot for the genome read mapped only plotting all families
ggplot((subset(TEfreq, repl== "0")) , aes(x = species, y = frac, fill=TEfam, order = -as.numeric(frac))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c(as.character(palette$col))) +
  scale_y_continuous(labels = percent_format()) +
  scale_x_discrete(limits=c("4_Tbi","4_Tte","3_Tce","3_Tms","1_Tps","1_Tdi","2_Tcm","2_Tsi","5_Tpa","5_Tge")) + 
  labs(x="species", y="fraction of reads mapped") +
  theme(axis.title = element_text(family = "Arial", color="#666666", face="bold", size=18), axis.text.x = element_text(family = "Arial", color="#666666", face="bold", size=14), axis.text.y = element_text(family = "Arial", color="#666666", face="bold", size=12)) +
  theme(panel.border = element_rect(linetype = "solid", colour = "grey", fill = NA), panel.grid.major = element_line(color = "grey", linetype = "dotted"), panel.grid.minor = element_line(colour = "grey", linetype = "dotted"), panel.background = element_blank(), axis.line = element_line(colour = "grey40")) +
  theme(legend.position = "right")

aggregate(TEfreq$frac, by=list(paste(TEfreq$species,TEfreq$repl)), FUN=sum)

#plot reseq one by one
ggplot(TEfreq, aes(x = (paste(species, "_", repl)), y = frac, fill=TEfam, order = -as.numeric(frac))) +
#ggplot(TEfreq, aes(x = species, y = frac, fill=TEfam, order = -as.numeric(frac))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c(as.character(palette$col))) +
  scale_y_continuous(labels = percent_format()) +
  #scale_x_discrete(limits=c("4_Tte_0","4_Tte_1","4_Tte_2","4_Tte_3","4_Tte_4","4_Tte_5","4_Tbi_0", "4_Tbi_1", "4_Tbi_2","4_Tbi_3","4_Tbi_4","4_Tbi_5")) + 
  labs(x="species", y="fraction of reads mapped") +
  theme(axis.title = element_text(family = "Arial", color="#666666", face="bold", size=18), axis.text.x = element_text(family = "Arial", color="#666666", face="bold", size=14), axis.text.y = element_text(family = "Arial", color="#666666", face="bold", size=12)) +
  theme(panel.border = element_rect(linetype = "solid", colour = "grey", fill = NA), panel.grid.major = element_line(color = "grey", linetype = "dotted"), panel.grid.minor = element_line(colour = "grey", linetype = "dotted"), panel.background = element_blank(), axis.line = element_line(colour = "grey40")) +
  theme(legend.position = "right") +
  theme(axis.text.x = element_text(angle = 90))


#add pair info
TEfreq$species2 <- ifelse(grepl("^4_Tbi$",TEfreq$species), "T. bartmani",
                          ifelse(grepl("4_Tte",TEfreq$species), "T. tahoe",
                                 ifelse(grepl("2_Tsi",TEfreq$species), "T. shepardi",
                                        ifelse(grepl("2_Tcm",TEfreq$species), "T. californicum",
                                               ifelse(grepl("3_Tms",TEfreq$species), "T. monikensis",
                                                      ifelse(grepl("3_Tce",TEfreq$species), "T. cristinae",
                                                             ifelse(grepl("1_Tdi",TEfreq$species), "T. douglasi",
                                                                    ifelse(grepl("1_Tps",TEfreq$species), "T. poppensis",
                                                                           ifelse(grepl("5_Tge",TEfreq$species), "T. genevievae",
                                               ifelse(grepl("5_Tpa",TEfreq$species), "T. podura", "Other"))))))))))
TEfreq$species2 <- as.factor(TEfreq$species2)

#plot reseq one by one
ggplot(TEfreq, aes(x = (paste(species2,repl)), y = frac, fill=TEfam, order = -as.numeric(frac))) +
  #ggplot(TEfreq, aes(x = species, y = frac, fill=TEfam, order = -as.numeric(frac))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c(as.character(palette$col))) +
  scale_y_continuous(labels = percent_format()) +
  scale_x_discrete(limits=c("T. bartmani 0","T. bartmani 1","T. bartmani 2","T. bartmani 3","T. bartmani 4","T. bartmani 5","T. tahoe 0","T. tahoe 1","T. tahoe 2","T. tahoe 3","T. tahoe 4","T. tahoe 5","T. shepardi 0","T. shepardi 3","T. shepardi 4","T. shepardi 5","T. californicum 0","T. californicum 1", "T. californicum 2","T. californicum 3","T. californicum 4","T. californicum 5","T. monikensis 0","T. monikensis 1","T. monikensis 2","T. monikensis 3","T. monikensis 4","T. monikensis 5","T. cristinae 0","T. cristinae 1","T. cristinae 2","T. cristinae 3","T. cristinae 4","T. cristinae 5","T. douglasi 0","T. douglasi 1","T. douglasi 2","T. douglasi 3","T. douglasi 4","T. douglasi 5","T. poppensis 0","T. poppensis 1","T. poppensis 2","T. poppensis 3","T. poppensis 4","T. poppensis 5","T. genevievae 0","T. genevievae 1","T. genevievae 2","T. genevievae 3","T. genevievae 4","T. genevievae 5","T. podura 0","T. podura 1","T. podura 2","T. podura 3","T. podura 4","T. podura 5")) + 
  labs(x="species", y="fraction of reads mapped") +
  theme(axis.title = element_text(family = "Arial", color="#666666", face="bold", size=18), axis.text.x = element_text(family = "Arial", color="#666666", face="italic", size=14), axis.text.y = element_text(family = "Arial", color="#666666", face="bold", size=12)) +
  theme(panel.border = element_rect(linetype = "solid", colour = "grey", fill = NA), panel.grid.major = element_line(color = "grey", linetype = "dotted"), panel.grid.minor = element_line(colour = "grey", linetype = "dotted"), panel.background = element_blank(), axis.line = element_line(colour = "grey40")) +
  theme(legend.position = "right") +
  theme(axis.text.x = element_text(angle = 90))




#Plotting of the total TE load including replicates (and 95% conf interval)
#Using median pairwise-dist age rank order

#MEAN
ggplot(TEfreq_agg, aes(x=species, y=frac))+ 
  stat_summary(geom="errorbar", fun.data=mean_cl_normal, width=0.4, fun.args=list(conf.int=0.95), color = "grey50") +
  stat_summary(geom="point", fun=mean, size=5, color=c("steelblue","orangered3","orangered3","steelblue","orangered3","steelblue","orangered3","steelblue","steelblue","orangered3")) +
  scale_y_continuous(labels = percent_format()) +
  scale_x_discrete(limits=c("4_Tbi","4_Tte","2_Tcm","2_Tsi","1_Tps","1_Tdi","3_Tce","3_Tms","5_Tpa","5_Tge")) +
  labs(x="species", y="fraction of reads mapped") +
  coord_cartesian(ylim=c(0,0.25)) +
  theme(axis.title = element_text(family = "Arial", color="#666666", face="bold", size=18), axis.text.x = element_text(family = "Arial", color="#666666", face="bold", size=14), axis.text.y = element_text(family = "Arial", color="#666666", face="bold", size=12)) +
  theme(panel.border = element_rect(linetype = "solid", colour = "grey", fill = NA), panel.grid.major = element_line(color = "grey", linetype = "dotted"), panel.grid.minor = element_line(colour = "grey", linetype = "dotted"), panel.background = element_blank(), axis.line = element_line(colour = "grey40")) +
  theme(legend.position = "right")

#MEDIAN
ggplot(TEfreq_agg, aes(x=species, y=frac))+ 
  stat_summary(geom="errorbar", fun.data=median_cl_boot, width=0.4, color = "grey50") +
  stat_summary(geom="point", fun=median, size=5, color=c("steelblue","orangered3","orangered3","steelblue","orangered3","steelblue","orangered3","steelblue","steelblue","orangered3")) +
  scale_y_continuous(labels = percent_format()) +
  scale_x_discrete(limits=c("4_Tbi","4_Tte","2_Tcm","2_Tsi","1_Tps","1_Tdi","3_Tce","3_Tms","5_Tpa","5_Tge")) +
  labs(x="species", y="fraction of reads mapped") +
  coord_cartesian(ylim=c(0,0.25)) +
  theme(axis.title = element_text(family = "Arial", color="#666666", face="bold", size=18), axis.text.x = element_text(family = "Arial", color="#666666", face="bold", size=14), axis.text.y = element_text(family = "Arial", color="#666666", face="bold", size=12)) +
  theme(panel.border = element_rect(linetype = "solid", colour = "grey", fill = NA), panel.grid.major = element_line(color = "grey", linetype = "dotted"), panel.grid.minor = element_line(colour = "grey", linetype = "dotted"), panel.background = element_blank(), axis.line = element_line(colour = "grey40")) +
  theme(legend.position = "right")



#Boxplot (mind the scale)
#ggplot(TEfreq_agg, aes(x = species, y = frac, fill=species, order = -as.numeric(frac))) +
#  stat_boxplot(geom ='errorbar') +
#  geom_boxplot() +
#  #geom_bar(stat = "identity") +
#  scale_fill_manual(values=c("steelblue","orangered3","orangered3","steelblue","orangered3","steelblue","orangered3","steelblue","steelblue","orangered3")) +
#  #scale_y_continuous(labels = percent_format()) +
#  scale_x_discrete(limits=c("4_Tbi","4_Tte","2_Tcm","2_Tsi","1_Tps","1_Tdi","3_Tce","3_Tms","5_Tpa","5_Tge")) +
#  labs(x="species", y="fraction of reads mapped") +
#  coord_cartesian(ylim=c(0.15,0.25)) +
#  theme(axis.title = element_text(family = "Arial", color="#666666", face="bold", size=18), axis.text.x = element_text(family = "Arial", color="#666666", face="bold", size=14), axis.text.y = element_text(family = "Arial", color="#666666", face="bold", size=12)) +
#  theme(panel.border = element_rect(linetype = "solid", colour = "grey", fill = NA), panel.grid.major = element_line(color = "grey", linetype = "dotted"), panel.grid.minor = element_line(colour = "grey", linetype = "dotted"), panel.background = element_blank(), axis.line = element_line(colour = "grey40")) +
#  theme(legend.position = "right")




#per TE family plot


#Example hwo to use facet_wrap
#ggplot(dat, aes(letters,value, label = letters)) + 
#  geom_bar(stat="identity") + 
#  facet_wrap(~grouping, scales="free")

#MEAN for each replicate split by TE superfamily and ordered by species pair
ggplot(TEfreq, aes(species, frac, label = species)) + 
  stat_summary(geom="errorbar", fun.data=mean_cl_normal, width=0.4, fun.args=list(conf.int=0.95), color = "grey50") +
  stat_summary(geom="point", fun=mean, size=2) +
  #geom_bar(stat="identity") + 
  facet_wrap(~TEfam, scales="free") +
  scale_y_continuous(labels = percent_format()) +
  scale_x_discrete(limits=c("4_Tbi","4_Tte","2_Tcm","2_Tsi","1_Tps","1_Tdi","3_Tce","3_Tms","5_Tpa","5_Tge")) +
  labs(x="species", y="fraction of reads mapped") +
  #coord_cartesian(ylim=c(0,0.055)) +
  theme(axis.title = element_text(family = "Arial", color="#666666", face="bold", size=12), axis.text.x = element_text(family = "Arial", color="#666666", face="bold", size=8), axis.text.y = element_text(family = "Arial", color="#666666", face="bold", size=8)) +
  theme(panel.border = element_rect(linetype = "solid", colour = "grey", fill = NA), panel.grid.major = element_line(color = "grey", linetype = "dotted"), panel.grid.minor = element_line(colour = "grey", linetype = "dotted"), panel.background = element_blank(), axis.line = element_line(colour = "grey40")) +
  theme(legend.position = "right")








#To combine super-families into Classes

#bin unclear TE superfamily assignment
TEfreq$TEfam <- sapply(TEfreq$TEfam, function(x){gsub(pattern = "DTX|DXX|RIX|RLX|RXX", replacement = "Unclear", x)}) #unclear TE supefamily assignment
#TEfreq$TEfam <- sapply(TEfreq$TEfam, function(x){gsub(pattern = "DXX|RXX", replacement = "unclear", x)}) #unclear order

TEfreq$TEfam <- as.factor(TEfreq$TEfam)
TEfreq$TEfam <- sapply(TEfreq$TEfam, function(x){gsub(pattern = "DTA|DTB|DTH|DTP|DTR|DTT|DTX", replacement = "TIR", x)})
TEfreq$TEfam <- sapply(TEfreq$TEfam, function(x){gsub(pattern = "RII|RIJ|RIL|RIR|RIT|RIX", replacement = "LINE", x)})
TEfreq$TEfam <- sapply(TEfreq$TEfam, function(x){gsub(pattern = "RLB|RLC|RLG|RLX", replacement = "LTR", x)})
TEfreq$TEfam <- sapply(TEfreq$TEfam, function(x){gsub(pattern = "RPP", replacement = "Penelope", x)})
TEfreq$TEfam <- sapply(TEfreq$TEfam, function(x){gsub(pattern = "DHH", replacement = "Helitron", x)})
TEfreq$TEfam <- sapply(TEfreq$TEfam, function(x){gsub(pattern = "DMM", replacement = "Maverick", x)})
TEfreq$TEfam <- sapply(TEfreq$TEfam, function(x){gsub(pattern = "RYX", replacement = "DIRS", x)})
TEfreq$TEfam <- sapply(TEfreq$TEfam, function(x){gsub(pattern = "RSX", replacement = "SINE", x)})
#TEfreq$TEfam <- sapply(TEfreq$TEfam, function(x){gsub(pattern = "RSX|XXX", replacement = "SINE", x)})

TEfreq$TEfam <- as.factor(TEfreq$TEfam)

str(TEfreq)


#sum count for the species and replicates levels (summing total TE load per replicate and species)
TEfreq_fam_agg <- aggregate(count~TEfam+species+repl, TEfreq, sum)
Mapped_fam_agg <- aggregate(mapped~TEfam+species+repl, TEfreq, mean)
TEfreq_fam_agg$mapped <- as.integer(Mapped_fam_agg$mapped)
TEfreq_fam_agg$frac <- TEfreq_fam_agg$count / TEfreq_fam_agg$mapped
#str(TEfreq_fam_agg)


#add repr mode info in new column
TEfreq_fam_agg$mode <- ifelse(grepl(match_asex,TEfreq_fam_agg$species), "asex",
                          ifelse(grepl(match_sex,TEfreq_fam_agg$species), "sex","Other"))
#add pair info
TEfreq_fam_agg$pair <- ifelse(grepl("Tdi|Tps",TEfreq_fam_agg$species), "Tdi-Tps",
                          ifelse(grepl("Tsi|Tcm",TEfreq_fam_agg$species), "Tsi-Tcm",
                                 ifelse(grepl("Tms|Tce",TEfreq_fam_agg$species), "Tms-Tce",
                                        ifelse(grepl("Tge|Tpa",TEfreq_fam_agg$species), "Tge-Tpa",
                                               ifelse(grepl("Tte|Tbi",TEfreq_fam_agg$species), "Tte-Tbi", "Other")))))

#modify data type
TEfreq_fam_agg$pair <- as.factor(TEfreq_fam_agg$pair)
TEfreq_fam_agg$mode <- as.factor(TEfreq_fam_agg$mode)
#TEfreq_fam_agg$frac <- as.numeric(TEfreq_fam_agg$frac)
TEfreq_fam_agg$repl <- as.factor(TEfreq_fam_agg$repl)

#change order
TEfreq_fam_agg$TEfam <- factor(TEfreq_fam_agg$TEfam, levels = c("TIR", "Helitron", "Maverick", "LTR", "DXX","DIRS", "Penelope", "LINE", "RXX","SINE","Unknown", "Unclear"))

#exclude some families
TEfreq_fam_agg <- subset(TEfreq_fam_agg, TEfam!= "RXX") #exclude unknown retro
TEfreq_fam_agg <- subset(TEfreq_fam_agg, TEfam!= "DXX") #exclude unknown retro


#MEAN for each replicate split by TE Class and ordered by species pair
ggplot(TEfreq_fam_agg, aes(species, frac, label = species)) + 
  stat_summary(geom="errorbar", fun.data=mean_cl_normal, width=0.4, fun.args=list(conf.int=0.95), color = "grey50") +
  stat_summary(geom="point", fun=mean, size=3) +
  #geom_bar(stat="identity") + 
  facet_wrap(~TEfam, scales="free") +
  scale_y_continuous(labels = percent_format()) +
  scale_x_discrete(limits=c("4_Tbi","4_Tte","2_Tcm","2_Tsi","1_Tps","1_Tdi","3_Tce","3_Tms","5_Tpa","5_Tge")) +
  labs(x="species", y="fraction of reads mapped") +
  #coord_cartesian(ylim=c(0,0.08)) +
  theme(axis.title = element_text(family = "Arial", color="#666666", face="bold", size=16), axis.text.x = element_text(family = "Arial", color="#666666", face="bold", size=12), axis.text.y = element_text(family = "Arial", color="#666666", face="bold", size=12)) +
  theme(panel.border = element_rect(linetype = "solid", colour = "grey", fill = NA), panel.grid.major = element_line(color = "grey", linetype = "dotted"), panel.grid.minor = element_line(colour = "grey", linetype = "dotted"), panel.background = element_blank(), axis.line = element_line(colour = "grey40")) +
  theme(legend.position = "right")



#TEclasses_samescale



#######
#STATS#
#######



#ANOVA (add pair in ANOVA plus interaction (~TE_type + pair * mode)) with 5000 bootstraps
#anovsF3(TEfreq$frac,TEfreq$TEfam,TEfreq$pair,TEfreq$mode)


#dev.off()
















#######################################################
#NOT IMPORTANT TESTS OF PLOTTING // REMOVED IN THE END#
########################################################

library(dplyr)

TEfreq_summary <- TEfreq_agg %>% # the names of the new data frame and the data frame to be summarised
  group_by(species) %>%   # the grouping variable
  summarise(mean_frac = mean(frac),  # calculates the mean of each group
            sd_frac = sd(frac), # calculates the standard deviation of each group
            n_frac = n(),  # calculates the sample size per group
            SE_frac = sd(frac)/sqrt(n())) # calculates the standard error of each group


ggplot(TEfreq_summary, aes(species, mean_frac)) + 
  geom_col() +  
  geom_errorbar(aes(ymin = mean_frac - sd_frac, ymax = mean_frac + sd_frac), width=0.2) +
  labs(y="frac ± s.d.", x = "Species") + theme_classic()



#following http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
#maybe for including the TEfams 
#TEfam_load <- summarySEwithin(TEfreq, measurevar="frac", betweenvars="species" ,withinvars="species",
#                         idvar="repl", na.rm=FALSE, conf.interval=.95)


TEfam_load_sum <- summarySEwithin(TEfreq_agg, measurevar="frac", withinvars="species", idvar="repl")

ggplot(TEfam_load_sum, aes(x=species, y=frac)) +
  geom_bar(position=position_dodge(.9), stat="identity", fill=c("steelblue","orangered3","orangered3","steelblue","orangered3","steelblue","orangered3","steelblue","steelblue","orangered3")) +
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=frac-ci, ymax=frac+ci)) +
  scale_fill_manual(values=c(as.character(palette$col))) +
  scale_y_continuous(labels = percent_format()) +
  scale_x_discrete(limits=c("4_Tbi","4_Tte","2_Tcm","2_Tsi","1_Tps","1_Tdi","3_Tce","3_Tms","5_Tpa","5_Tge")) + 
  labs(x="species", y="fraction of reads mapped") +
  theme(axis.title = element_text(family = "Arial", color="#666666", face="bold", size=18), axis.text.x = element_text(family = "Arial", color="#666666", face="bold", size=14), axis.text.y = element_text(family = "Arial", color="#666666", face="bold", size=12)) +
  theme(panel.border = element_rect(linetype = "solid", colour = "grey", fill = NA), panel.grid.major = element_line(color = "grey", linetype = "dotted"), panel.grid.minor = element_line(colour = "grey", linetype = "dotted"), panel.background = element_blank(), axis.line = element_line(colour = "grey40")) +
  theme(legend.position = "right")


#multiple levels (TE families)
datac <- summarySEwithin(TEfreq, measurevar="frac", withinvars=c("species","TEfam"), idvar="repl")


ggplot(datac, aes(x=species, y=frac, label=species)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=frac-sd, ymax=frac+sd)) +
  facet_wrap(~TEfam, scales="free") +
  scale_fill_manual(values=c(as.character(palette$col))) +
  scale_y_continuous(labels = percent_format()) +
  scale_x_discrete(limits=c("4_Tbi","4_Tte","2_Tcm","2_Tsi","1_Tps","1_Tdi","3_Tce","3_Tms","5_Tpa","5_Tge")) + 
  labs(x="species", y="fraction of reads mapped") +
  theme(axis.title = element_text(family = "Arial", color="#666666", face="bold", size=18), axis.text.x = element_text(family = "Arial", color="#666666", face="bold", size=14), axis.text.y = element_text(family = "Arial", color="#666666", face="bold", size=12)) +
  theme(panel.border = element_rect(linetype = "solid", colour = "grey", fill = NA), panel.grid.major = element_line(color = "grey", linetype = "dotted"), panel.grid.minor = element_line(colour = "grey", linetype = "dotted"), panel.background = element_blank(), axis.line = element_line(colour = "grey40")) +
  theme(legend.position = "right")



##################
#HELPER FUNCTIONS#
##################


## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}




## Norms the data within specified groups in a data frame; it normalizes each
## subject (identified by idvar) so that they have the same mean, within each group
## specified by betweenvars.
##   data: a data frame.
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   na.rm: a boolean that indicates whether to ignore NA's
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
  library(plyr)
  
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                         .fun = function(xx, col, na.rm) {
                           c(subjMean = mean(xx[,col], na.rm=na.rm))
                         },
                         measurevar,
                         na.rm
  )
  
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}




## Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, un-normed mean, normed mean (with same between-group mean),
##   standard deviation, standard error of the mean, and confidence interval.
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                       FUN=is.factor, FUN.VALUE=logical(1))
  
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  
  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  
  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                  FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}











