### positive_sel_plot_poly.R

library(ggplot2)
library(cowplot)
library(hash)
library(stringr)
library(lme4)
library(car)
library(glmmTMB)

library("arm")
library("lme4")
library("MuMIn")

print (sessionInfo())
# 
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Big Sur 10.16
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] glmmTMB_1.0.2.1 MuMIn_1.43.17   arm_1.11-2      MASS_7.3-53     car_3.0-10     
# [6] carData_3.0-4   lme4_1.1-26     Matrix_1.3-2    stringr_1.4.0   hash_2.2.6.1   
# [11] cowplot_1.1.1   ggplot2_3.3.3  
# 
# loaded via a namespace (and not attached):
#   [1] splines_4.0.3       Formula_1.2-4       statmod_1.4.35      stats4_4.0.3       
# [5] latticeExtra_0.6-29 cellranger_1.1.0    pillar_1.4.7        backports_1.2.1    
# [9] lattice_0.20-41     glue_1.4.2          digest_0.6.27       RColorBrewer_1.1-2 
# [13] checkmate_2.0.0     minqa_1.2.4         sandwich_3.0-0      colorspace_2.0-0   
# [17] htmltools_0.5.1.1   pkgconfig_2.0.3     haven_2.3.1         purrr_0.3.4        
# [21] xtable_1.8-4        mvtnorm_1.1-1       scales_1.1.1        jpeg_0.1-8.1       
# [25] openxlsx_4.2.3      rio_0.5.16          emmeans_1.6.1       htmlTable_2.1.0    
# [29] tibble_3.0.6        generics_0.1.0      farver_2.0.3        ellipsis_0.3.1     
# [33] TH.data_1.0-10      withr_2.4.1         nnet_7.3-15         TMB_1.7.20         
# [37] survival_3.2-7      magrittr_2.0.1      crayon_1.4.0        readxl_1.3.1       
# [41] estimability_1.3    nlme_3.1-151        forcats_0.5.1       foreign_0.8-81     
# [45] tools_4.0.3         data.table_1.13.6   hms_1.0.0           multcomp_1.4-15    
# [49] lifecycle_0.2.0     munsell_0.5.0       cluster_2.1.0       zip_2.1.1          
# [53] compiler_4.0.3      rlang_0.4.10        grid_4.0.3          nloptr_1.2.2.2     
# [57] rstudioapi_0.13     htmlwidgets_1.5.3   base64enc_0.1-3     labeling_0.4.2     
# [61] boot_1.3-26         codetools_0.2-18    gtable_0.3.0        abind_1.4-5        
# [65] curl_4.3            R6_2.5.0            zoo_1.8-8           gridExtra_2.3      
# [69] knitr_1.31          dplyr_1.0.3         Hmisc_4.4-2         stringi_1.5.3      
# [73] Rcpp_1.0.6          vctrs_0.3.6         rpart_4.1-15        png_0.1-7          
# [77] tidyselect_1.1.0    xfun_0.20           coda_0.19-4 

#############################################################################
###

## data
poly_dat <- read.csv("pos_sel_data/allsp_poly_or_fixed.csv") ## made with pos_poly.py
head(poly_dat)

### get genes with non-syn polymorphic sites
poly_HOGS_a <- unique(subset(poly_dat, poly_dat$poly_aa == "poly")$HOG)

### also remove sister sp HOG of genes with non-syn polymorphic sites to keep number of HOGs in sex-asex pairs the same
poly_HOGS <- c()
for(i in poly_HOGS_a){
  print(i)
  curr_sp  <- strsplit(i, "__")[[1]][1]
  curr_HOG <- strsplit(i, "__")[[1]][2]
  print(curr_sp)
  new_sp = ""
  if(curr_sp == "Tbi"){new_sp = "Tte"}
  if(curr_sp == "Tte"){new_sp = "Tbi"}
  
  if(curr_sp == "Tce"){new_sp = "Tms"}
  if(curr_sp == "Tms"){new_sp = "Tce"}
  
  if(curr_sp == "Tcm"){new_sp = "Tsi"}
  if(curr_sp == "Tsi"){new_sp = "Tcm"}
  
  if(curr_sp == "Tpa"){new_sp = "Tge"}
  if(curr_sp == "Tge"){new_sp = "Tpa"}
  
  if(curr_sp == "Tps"){new_sp = "Tdi"}
  if(curr_sp == "Tdi"){new_sp = "Tps"}
  
  print(new_sp)
  
  poly_HOGS <- c(poly_HOGS, i, paste(new_sp, curr_HOG, sep = "__"))
}

### remove poly genes from selection analyses
dat1      <- read.table("pos_sel_data/timema_543_branches_with-ncat-codon-rate_sites_with_h0.tsv", sep = "\t", header = T)
dat1$gene <- as.character(dat1$gene )
dat1$sp_gene <- paste(dat1$branch_name, dat1$gene, sep = "__" )
length(dat1[,1])
dat1 = subset(dat1, ! dat1$sp_gene %in% poly_HOGS)
length(dat1[,1])
head(dat1)


dat1$branch_name <- as.character(dat1$branch_name)
dat1$branch_name <- 
ifelse(dat1$branch_name == "Northern_Clade", "Northern", 
ifelse(dat1$branch_name == "Santa_Barbara_Clade", "Santa Barbara", 
ifelse(dat1$branch_name == "Southern_Clade", "Southern", 
dat1$branch_name)))

dat1$sp_pair <-
ifelse(dat1$branch_name == "Tbi", "Tbi-Tte",
ifelse(dat1$branch_name == "Tce", "Tce-Tms",
ifelse(dat1$branch_name == "Tcm", "Tcm-Tsi",
ifelse(dat1$branch_name == "Tpa", "Tpa-Tge",
ifelse(dat1$branch_name == "Tps", "Tps-Tdi",
ifelse(dat1$branch_name == "Tte", "Tbi-Tte", 
ifelse(dat1$branch_name == "Tms", "Tce-Tms",
ifelse(dat1$branch_name == "Tsi", "Tcm-Tsi",
ifelse(dat1$branch_name == "Tge", "Tpa-Tge",
ifelse(dat1$branch_name == "Tdi", "Tps-Tdi", 
NA))))))))))

dat1$rep_mode <-
ifelse(dat1$branch_name == "Tbi", "sex",
ifelse(dat1$branch_name == "Tce", "sex",
ifelse(dat1$branch_name == "Tcm", "sex",
ifelse(dat1$branch_name == "Tpa", "sex",
ifelse(dat1$branch_name == "Tps", "sex",
ifelse(dat1$branch_name == "Tte", "asex", 
ifelse(dat1$branch_name == "Tms", "asex", 
ifelse(dat1$branch_name == "Tsi", "asex",
ifelse(dat1$branch_name == "Tge", "asex",
ifelse(dat1$branch_name == "Tdi", "asex", 
ifelse(dat1$branch_name == "Tps/Tdi", "sex_asex", 
ifelse(dat1$branch_name == "Tpa/Tge", "sex_asex", 
ifelse(dat1$branch_name == "Tcm/Tsi", "sex_asex", 
ifelse(dat1$branch_name == "Tbi/Tte", "sex_asex", 
ifelse(dat1$branch_name == "Santa Barbara", "sex_asex", 
ifelse(dat1$branch_name == "Northern", "clade",
ifelse(dat1$branch_name == "Southern", "clade",
NA)))))))))))))))))

head(dat1)

#######################################################################
### how many branches show +ve sel by sp

q_val_threshold <- 0.01
dat1_selected = subset(dat1, dat1$qvalue < q_val_threshold)
length(dat1_selected[,1])
head(dat1_selected )
min(dat1_selected$lrt)

N_sel_branches <- as.data.frame(table(dat1_selected$branch_name))
colnames(N_sel_branches) <- c("branch", "N")

N_sel_branches$rep_mode <-
ifelse(N_sel_branches$branch == "Tbi", "sex",
ifelse(N_sel_branches$branch == "Tce", "sex",
ifelse(N_sel_branches$branch == "Tcm", "sex",
ifelse(N_sel_branches$branch == "Tpa", "sex",
ifelse(N_sel_branches$branch == "Tps", "sex",
ifelse(N_sel_branches$branch == "Tte", "asex", 
ifelse(N_sel_branches$branch == "Tms", "asex", 
ifelse(N_sel_branches$branch == "Tsi", "asex",
ifelse(N_sel_branches$branch == "Tge", "asex",
ifelse(N_sel_branches$branch == "Tdi", "asex", 
ifelse(N_sel_branches$branch == "Tps/Tdi", "sex_asex", 
ifelse(N_sel_branches$branch == "Tpa/Tge", "sex_asex", 
ifelse(N_sel_branches$branch == "Tcm/Tsi", "sex_asex", 
ifelse(N_sel_branches$branch == "Tbi/Tte", "sex_asex", 
ifelse(N_sel_branches$branch == "Santa Barbara", "sex_asex", 
ifelse(N_sel_branches$branch == "Northern", "clade",
ifelse(N_sel_branches$branch == "Southern", "clade",
NA)))))))))))))))))

N_sel_branches_without_NandS    <- subset(N_sel_branches, N_sel_branches$rep_mode != "clade")
N_sel_branches_without_termonly <- subset(N_sel_branches_without_NandS, N_sel_branches_without_NandS$rep_mode != "sex_asex")
N_sel_branches_without_termonly$group <- c("Tbi-Tte", "Tce-Tms", "Tcm-Tsi", "Tps-Tdi", "Tpa-Tge", "Tce-Tms", "Tpa-Tge", "Tps-Tdi", "Tcm-Tsi", "Tbi-Tte")

N_sel_branches_without_termonly <- subset(N_sel_branches_without_NandS, N_sel_branches_without_NandS$rep_mode != "sex_asex")
N_sel_branches_without_termonly$group <- c("Tbi-Tte", "Tce-Tms", "Tcm-Tsi", "Tps-Tdi", "Tpa-Tge", "Tce-Tms", "Tpa-Tge", "Tps-Tdi", "Tcm-Tsi", "Tbi-Tte")
N_sel_branches_without_termonly$rep_mode_ord  = ordered(N_sel_branches_without_termonly$rep_mode, levels = c("asex", "sex"))
N_sel_branches_without_termonly$group_ord = ordered(N_sel_branches_without_termonly$group, levels = c("Tbi-Tte", "Tcm-Tsi", "Tce-Tms", "Tps-Tdi", "Tpa-Tge"))

max_y = max(N_sel_branches$N * 1.05)

P1b <- ggplot(N_sel_branches_without_termonly, aes(x = factor(group_ord), y = N, fill = rep_mode_ord)) + 
	geom_col(width = 0.5, colour="black", position=position_dodge(width=0.6)) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
	scale_fill_manual(values = c("#92C5DECD", "#D6604DED")) + 
	xlab ("Species pair") + 
	ylab ("Number of positive selected genes")  + 
	ggtitle(paste("543sp pair 1-to-1 orths, BSG", "qval thresh = ",q_val_threshold ))  + ylim(0,max_y )


pdf(paste("543sp_N_pos_genes_int_branches_P1b_qval_", q_val_threshold ,"_poly_aa_filt.pdf", sep = ""), width = 6, height = 8)
P1b  
dev.off()
getwd() ## where has my plot gone....



# png(filename = paste("N_pos_genes_int_branches_P3a_", q_val_threshold ,".png", sep = ""), width = 8, height = 8, units = "in", bg = "white", res = 300)
# P3a 
# dev.off()
# getwd() ## where has my plot gone....




###################################################################################################
### sig diff

dat1_a    <- subset(dat1, dat1$rep_mode != "clade")
dat1_term <- subset(dat1_a, dat1_a$rep_mode != "sex_asex")
dat1_term$pos_sel_bi <- ifelse(dat1_term$qvalue < q_val_threshold, 1, 0)
length(dat1_term[,1])
head(dat1_term)
dat1_term$rep_mode <- as.factor(dat1_term$rep_mode)
dat1_term$gene <- as.factor(dat1_term$gene)
str(dat1_term)
length(na.omit(dat1_term)[,1])
length(dat1_term[,1])


### rep mode inter sp_pair

### test interaction	
mix_5 = glmer(pos_sel_bi ~ sp_pair * rep_mode + (1|gene), family = "binomial", control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = dat1_term) #### Runs
Anova(mix_5, type = 3 )  ### WALD

### test main effects	
mix_5a = glmer(pos_sel_bi ~ sp_pair + rep_mode + (1|gene), family = "binomial", control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e9)), data = dat1_term) #### fails for some values of q_val_threshold. in which case use TMB
Anova(mix_5a, type = 3 )  ### WALD 


### test interaction	TMB
mix_5_TMB = glmmTMB(pos_sel_bi ~ sp_pair * rep_mode + (1|gene), family = "binomial", data = dat1_term) #### Runs
Anova(mix_5_TMB, type = 3 )  ### WALD

mix_5a_TMB = glmmTMB(pos_sel_bi ~ sp_pair + rep_mode + (1|gene), family = "binomial", data = dat1_term) #### Runs
Anova(mix_5a_TMB, type = 3 )  ### WALD













