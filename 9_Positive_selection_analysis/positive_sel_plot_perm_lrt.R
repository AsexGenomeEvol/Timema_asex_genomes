### +ve sel

library(ggplot2)
library(cowplot)
library(hash)
library(stringr)
library(car)
library(MASS)
library(fitdistrplus)

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
 # [1] fitdistrplus_1.0-14 npsurv_0.4-0        lsei_1.2-0          survival_3.1-7      MASS_7.3-51.4       car_3.0-3           carData_3.0-2       stringr_1.4.0       hash_2.2.6.1        cowplot_1.0.0       ggplot2_3.3.2      

# loaded via a namespace (and not attached):
 # [1] zip_2.0.3         Rcpp_1.0.2        pillar_1.4.2      compiler_3.5.1    cellranger_1.1.0  forcats_0.4.0     tools_3.5.1       lattice_0.20-38   lifecycle_0.2.0   tibble_2.1.3      gtable_0.3.0      pkgconfig_2.0.2  
# [13] rlang_0.4.8       Matrix_1.2-17     openxlsx_4.1.0.1  curl_4.0          haven_2.1.1       rio_0.5.16        withr_2.1.2       dplyr_1.0.2       generics_0.0.2    vctrs_0.3.4       hms_0.5.1         grid_3.5.1       
# [25] tidyselect_1.1.0  glue_1.4.2        data.table_1.12.2 R6_2.4.0          readxl_1.3.1      foreign_0.8-72    purrr_0.3.2       magrittr_1.5      splines_3.5.1     scales_1.0.0      abind_1.4-5       colorspace_1.4-1 
# [37] stringi_1.4.3     munsell_0.5.0     crayon_1.3.4     


## data
dat1      <- read.table("pos_sel_data/timema_543_branches_with-ncat-codon-rate_sites_with_h0.tsv", sep = "\t", header = T)
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

###################################################################################################
### sig diff

dat1_a    <- subset(dat1, dat1$rep_mode != "clade")
dat1_term <- subset(dat1_a, dat1_a$rep_mode != "sex_asex")
dat1_term$rep_mode <- as.factor(dat1_term$rep_mode)
dat1_term$gene <- as.factor(dat1_term$gene)



#####################################
### get test stats from real data

### using 2 dists to check it is robust
m7c_real = glm(dat1_term$lrt ~ dat1_term$sp_pair + dat1_term$rep_mode)
m7d_real = glm(dat1_term$lrt ~ dat1_term$sp_pair + dat1_term$rep_mode, family = quasipoisson(link = "log"))

m7c_real_sp_pair_LR <- Anova(m7c_real, type = 3)$LR[1]
m7d_real_sp_pair_LR <- Anova(m7d_real, type = 3)$LR[1]

m7c_real_rep_mode_LR <- Anova(m7c_real, type = 3)$LR[2]
m7d_real_rep_mode_LR <- Anova(m7d_real, type = 3)$LR[2]


########## randomise rep mode

rand_rep_mode <- function(df){
	pos <- c("sex", "asex")
	rand_rep <- c()
	for (i in seq(1,length(df[,1]) / 2)){
		rand_rep_i <- sample(pos, replace = F)
		rand_rep <- c(rand_rep, rand_rep_i) 	
	}
	df$rand_rep <- rand_rep
	m7c = glm(df$lrt ~ df$sp_pair + df$rand_rep)
	m7c_out = Anova(m7c, type = 3)
	m7d = glm(df$lrt ~ df$sp_pair + df$rand_rep, family = quasipoisson(link = "log"))
	m7d_out = Anova(m7d, type = 3)
	
	## LR sp pair, LR rep, P sp pair, P rep mode
	m7c_out_v <- c(m7c_out$LR[1], m7c_out$LR[2], m7c_out$P[1], m7c_out$P[2])
	m7d_out_v <- c(m7d_out$LR[1], m7d_out$LR[2], m7d_out$P[1], m7d_out$P[2])
	output <- list("m7c_out_v" = m7c_out_v, "m7d_out_v" = m7d_out_v )
	
	return(output)	
}

#### run for x times

run_N = 1000 ### number of randomisations. this takes some time to run.
set.seed(42)

rand_rep_df_m7c <- c()
rand_rep_df_m7d <- c()
for (i in seq(1:run_N)){
	print(i)
	test_i <- rand_rep_mode(dat1_term)
	rand_rep_df_m7c <- rbind(rand_rep_df_m7c, test_i$m7c_out_v)
	rand_rep_df_m7d <- rbind(rand_rep_df_m7d, test_i$m7d_out_v)
}


colnames(rand_rep_df_m7c) <- c("LR_sp","LR_rep_mode", "P_sp","P_rep_mode")
rand_rep_df_m7c <- as.data.frame(rand_rep_df_m7c)

colnames(rand_rep_df_m7d) <- c("LR_sp","LR_rep_mode", "P_sp","P_rep_mode")
rand_rep_df_m7d <- as.data.frame(rand_rep_df_m7d)

get_pval = function(rand_df,calc_TS){
	N_rand_larger <- length(subset(rand_df, rand_df$LR_rep_mode > calc_TS)[,1])
	print(N_rand_larger)
	
	adj_pval = 10000
	if(N_rand_larger == 0){
		adj_pval = 0
	}
	else{
		adj_pval = N_rand_larger / length(rand_df$LR_rep_mode)
	}
	
	print(adj_pval)
}


get_pval(rand_rep_df_m7c,  m7c_real_sp_pair_LR ) # 0.011
get_pval(rand_rep_df_m7d,  m7d_real_sp_pair_LR ) # 0.011

 
