### +ve sel

library(ggplot2)
library(cowplot)
library(hash)
library(stringr)
library(lme4)
library(car)

print (sessionInfo())

# R version 3.5.1 (2018-07-02)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS  10.15.7

# Matrix products: default
# BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
# [1] car_3.0-3     carData_3.0-2 lme4_1.1-21   Matrix_1.2-17 stringr_1.4.0 hash_2.2.6.1  cowplot_1.0.0 ggplot2_3.3.2

# loaded via a namespace (and not attached):
 # [1] zip_2.0.3         Rcpp_1.0.2        cellranger_1.1.0  pillar_1.4.2      compiler_3.5.1    nloptr_1.2.1      forcats_0.4.0     tools_3.5.1       boot_1.3-23       lifecycle_0.2.0   tibble_2.1.3      gtable_0.3.0     
# [13] nlme_3.1-141      lattice_0.20-38   pkgconfig_2.0.2   rlang_0.4.8       openxlsx_4.1.0.1  curl_4.0          haven_2.1.1       rio_0.5.16        withr_2.1.2       dplyr_1.0.2       hms_0.5.1         generics_0.0.2   
# [25] vctrs_0.3.4       grid_3.5.1        tidyselect_1.1.0  glue_1.4.2        data.table_1.12.2 R6_2.4.0          readxl_1.3.1      foreign_0.8-72    minqa_1.2.4       purrr_0.3.2       magrittr_1.5      scales_1.0.0     
# [37] splines_3.5.1     MASS_7.3-51.4     abind_1.4-5       colorspace_1.4-1  stringi_1.4.3     munsell_0.5.0     crayon_1.3.4  

## data
dat1      <- read.table("data/timema_543_branches_with-ncat-codon-rate_sites_with_h0.tsv", sep = "\t", header = T)
dat1$gene <- as.character(dat1$gene )
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

q_val_threshold <- 0.05 ## set to desired threshold
dat1_selected = subset(dat1, dat1$qvalue < q_val_threshold)

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

pdf(paste("543sp_N_pos_genes_int_branches_P1b_qval_", q_val_threshold ,".pdf", sep = ""), width = 6, height = 8)
P1b  
dev.off()
getwd() ## where has my plot gone....



###################################################################################################
### sig diff


dat1_a    <- subset(dat1, dat1$rep_mode != "clade")
dat1_term <- subset(dat1_a, dat1_a$rep_mode != "sex_asex")

dat1_term$pos_sel_bi <- ifelse(dat1_term$qvalue < q_val_threshold, 1, 0)
head(dat1_term)
dat1_term$rep_mode <- as.factor(dat1_term$rep_mode)
dat1_term$gene <- as.factor(dat1_term$gene)
str(dat1_term)


attach(dat1_term)


### test interaction	
mix_5 = glmer(pos_sel_bi ~ sp_pair * rep_mode + (1|gene), family = "binomial", control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7))) 
Anova(mix_5, type = 3 )  ### WALD

### test main effects	
mix_5a = glmer(pos_sel_bi ~ sp_pair + rep_mode + (1|gene), family = "binomial", control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7))) 
Anova(mix_5a, type = 3 )  ### WALD


