### +ve sel

library(ggplot2)
library(cowplot)
library(hash)
library(stringr)
library("SuperExactTest")
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
# [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
# [1] SuperExactTest_1.0.7 stringr_1.4.0        hash_2.2.6.1         cowplot_1.0.0        ggplot2_3.3.2       

# loaded via a namespace (and not attached):
 # [1] Rcpp_1.0.2       withr_2.1.2      crayon_1.3.4     dplyr_1.0.2      R6_2.4.0         lifecycle_0.2.0  gtable_0.3.0     magrittr_1.5     scales_1.0.0     pillar_1.4.2     stringi_1.4.3    rlang_0.4.8      vctrs_0.3.4     
# [14] generics_0.0.2   tools_3.5.1      glue_1.4.2       purrr_0.3.2      munsell_0.5.0    compiler_3.5.1   pkgconfig_2.0.2  colorspace_1.4-1 tidyselect_1.1.0 tibble_2.1.3    

Intersect <- function (x) {  
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}


#####################
### overlaps between positively selected genes - sep as different N genes in each comp


dat1      <- read.table("pos_sel_data/timema_543_branches_with-ncat-codon-rate_sites_with_h0.tsv", sep = "\t", header = T)
dat1$gene <- as.character(dat1$gene )
head(dat1)


calc_overlap = function(df, sp_A, sp_B, qval){
	df_spA = subset(dat1, df$branch_name == sp_A)
	df_spB = subset(dat1, df$branch_name == sp_B)

	df_spA_sig = subset(df_spA, df_spA$qval < qval)
	df_spB_sig = subset(df_spB, df_spB$qval < qval)

	
	
	all_gene_A = df_spA$gene
	all_gene_B = df_spB$gene

	sig_gene_A = df_spA_sig$gene
	sig_gene_B = df_spB_sig$gene	
	
	print(Intersect(list(sig_gene_A, sig_gene_B)))
	
	all_gene_list = list(sp_A = all_gene_A, sp_B = all_gene_B)
	all_gene_fit_overlap  = length(MSET(all_gene_list, n=length(df[,1]), lower.tail=FALSE)$intersects)
	
	sig_gene_list  = list(sp_A = sig_gene_A, sp_B = sig_gene_B)
	all_gene_fit   = MSET(sig_gene_list, n=all_gene_fit_overlap, lower.tail=FALSE)	
	# print(all_gene_fit_overlap)
	
	# print(sig_gene_A)
	# print(all_gene_fit)
	
	
	# expected 
	length.gene.sets_sig =sapply(sig_gene_list,length)
	#print(length.gene.sets_sig)
	num.expcted.overlap=all_gene_fit_overlap*do.call(prod,as.list(length.gene.sets_sig/all_gene_fit_overlap))
	
	out_v <- c(length(all_gene_fit$intersects), num.expcted.overlap, all_gene_fit_overlap, all_gene_fit$FE, all_gene_fit$p)
	out_df <- as.data.frame(rbind(out_v))
	colnames(out_df) <- c("N_overlap", "N_exp", "Total_genes", "FE", "p")
	rownames(out_df) <- paste(sp_A, sp_B, sep = "_")
	return(out_df)
	
}


all_combo <- expand.grid(a = c("Tbi", "Tte", "Tce", "Tms","Tcm", "Tsi", "Tpa", "Tge", "Tps", "Tdi"), b = c("Tbi", "Tte", "Tce", "Tms","Tcm", "Tsi", "Tpa", "Tge", "Tps", "Tdi"), unique = TRUE)
all_combo$a <- as.character(all_combo$a)
all_combo$b <- as.character(all_combo$b)


###### get uniq combos

all_combo_list_s <- list()
for(i in seq(1: length(all_combo[,1]))){
	sp_sort = sort(c(all_combo$a[i], all_combo$b[i]))
	sp_str = paste(sp_sort[1], "_", sp_sort[2], sep = "")
	
	print(sp_sort)	
	all_combo_list_s <- c(all_combo_list_s, sp_str)
}

all_combo_list_s_uniq <- unique(all_combo_list_s)

#### run for all uniq combos
all_combo_overlap <- c()
for(i in seq(1: length(all_combo_list_s_uniq))){
	rec = as.character(all_combo_list_s_uniq[i])
	rec_s <- strsplit(rec, "_")
	print(rec)
	print(rec_s[[1]][1])
	print(rec_s[[1]][2])	
	over_i <- calc_overlap(dat1, rec_s[[1]][1], rec_s[[1]][2], 0.05)
	all_combo_overlap = rbind(all_combo_overlap, over_i)
}


#################
## remove species compared to itself (e.g Tbi Tbi)

all_combo_overlap$sp_a <- str_split_fixed(as.character(rownames(all_combo_overlap)), "_", 2)[,1]
all_combo_overlap$sp_b <- str_split_fixed(as.character(rownames(all_combo_overlap)), "_", 2)[,2]
all_combo_overlap$sp_same <- ifelse(all_combo_overlap$sp_a == all_combo_overlap$sp_b, 0, 1)

all_combo_overlap_want <- subset(all_combo_overlap, all_combo_overlap$sp_same == 1)

length(all_combo_overlap_want[,1]) ## 45

### FDR p-vals

all_combo_overlap_want$FDR <- p.adjust(all_combo_overlap_want$p, method = "fdr")
min(all_combo_overlap_want$FDR) ## 0.4940856


####### output 

write.csv(all_combo_overlap_want, "all_combo_overlap_want.csv", row.names = F, quote = F)



