### topGO
# install
# source("http://bioconductor.org/biocLite.R") 
# biocLite() 
# source("http://bioconductor.org/biocLite.R")   
# biocLite("topGO")
# biocLite("ALL")
# biocLite("affyLib")

library(topGO)
library(ALL)
library("VennDiagram")
library(gridExtra)
library(grid)
library(ggplot2)
library("SuperExactTest")
library(cowplot)
require(dplyr)
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
 # [1] grid      stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
 # [1] dplyr_1.0.2          cowplot_1.0.0        SuperExactTest_1.0.7 ggplot2_3.3.2        gridExtra_2.3        VennDiagram_1.6.20   futile.logger_1.4.3  ALL_1.24.0           topGO_2.34.0         SparseM_1.77         GO.db_3.7.0         
# [12] AnnotationDbi_1.44.0 IRanges_2.16.0       S4Vectors_0.20.1     Biobase_2.42.0       graph_1.60.0         BiocGenerics_0.28.0 

# loaded via a namespace (and not attached):
 # [1] Rcpp_1.0.2           compiler_3.5.1       pillar_1.4.2         formatR_1.7          futile.options_1.0.1 digest_0.6.20        bit_1.1-14           lifecycle_0.2.0      RSQLite_2.1.2        memoise_1.1.0        tibble_2.1.3        
# [12] gtable_0.3.0         lattice_0.20-38      pkgconfig_2.0.2      rlang_0.4.8          DBI_1.0.0            withr_2.1.2          generics_0.0.2       vctrs_0.3.4          tidyselect_1.1.0     bit64_0.9-7          glue_1.4.2          
# [23] R6_2.4.0             purrr_0.3.2          magrittr_1.5         lambda.r_1.2.3       blob_1.2.0           scales_1.0.0         matrixStats_0.54.0   colorspace_1.4-1     munsell_0.5.0        crayon_1.3.4          

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

Union <- function (x) {  
  # Multiple set version of union
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    union(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    union(x[[1]], Union(x[-1]))
  }
}

Setdiff <- function (x, y) {
  # Remove the union of the y's from the common x's. 
  # x and y are lists of characters.
  xx <- Intersect(x)
  yy <- Union(y)
  setdiff(xx, yy)
}

#### load annotation

setwd("data/GO_terms")

# ## Arth annotated ## not used
# geneID2GO_Tbi_Arth <- readMappings(file = "543sp_Arth_Tbi_forTopGO.txt")
# geneID2GO_Tte_Arth <- readMappings(file = "543sp_Arth_Tte_forTopGO.txt")
# geneID2GO_Tce_Arth <- readMappings(file = "543sp_Arth_Tce_forTopGO.txt")
# geneID2GO_Tms_Arth <- readMappings(file = "543sp_Arth_Tms_forTopGO.txt")
# geneID2GO_Tcm_Arth <- readMappings(file = "543sp_Arth_Tcm_forTopGO.txt")
# geneID2GO_Tsi_Arth <- readMappings(file = "543sp_Arth_Tsi_forTopGO.txt")
# geneID2GO_Tpa_Arth <- readMappings(file = "543sp_Arth_Tpa_forTopGO.txt")
# geneID2GO_Tge_Arth <- readMappings(file = "543sp_Arth_Tge_forTopGO.txt")
# geneID2GO_Tps_Arth <- readMappings(file = "543sp_Arth_Tps_forTopGO.txt")
# geneID2GO_Tdi_Arth <- readMappings(file = "543sp_Arth_Tdi_forTopGO.txt")


## Droso annotated
geneID2GO_Tbi_Droso <- readMappings(file = "543sp_Droso_Tbi_forTopGO.txt")
geneID2GO_Tte_Droso <- readMappings(file = "543sp_Droso_Tte_forTopGO.txt")
geneID2GO_Tce_Droso <- readMappings(file = "543sp_Droso_Tce_forTopGO.txt")
geneID2GO_Tms_Droso <- readMappings(file = "543sp_Droso_Tms_forTopGO.txt")
geneID2GO_Tcm_Droso <- readMappings(file = "543sp_Droso_Tcm_forTopGO.txt")
geneID2GO_Tsi_Droso <- readMappings(file = "543sp_Droso_Tsi_forTopGO.txt")
geneID2GO_Tpa_Droso <- readMappings(file = "543sp_Droso_Tpa_forTopGO.txt")
geneID2GO_Tge_Droso <- readMappings(file = "543sp_Droso_Tge_forTopGO.txt")
geneID2GO_Tps_Droso <- readMappings(file = "543sp_Droso_Tps_forTopGO.txt")
geneID2GO_Tdi_Droso <- readMappings(file = "543sp_Droso_Tdi_forTopGO.txt")



###############################################################################################################################################
#### read in tables with genename and qval



make_named_numeric_vector <- function(list_file_name){
	full_list <- as.list(read.table(list_file_name))
	full_list_GL <- full_list$V2
	names(full_list_GL) <- full_list$V1
	return(full_list_GL)
}

Tbi_qval_GL <- make_named_numeric_vector("543sp_Arth_Tbi_qval.txt")
Tte_qval_GL <- make_named_numeric_vector("543sp_Arth_Tte_qval.txt")
Tce_qval_GL <- make_named_numeric_vector("543sp_Arth_Tce_qval.txt")
Tms_qval_GL <- make_named_numeric_vector("543sp_Arth_Tms_qval.txt")
Tcm_qval_GL <- make_named_numeric_vector("543sp_Arth_Tcm_qval.txt")
Tsi_qval_GL <- make_named_numeric_vector("543sp_Arth_Tsi_qval.txt")
Tpa_qval_GL <- make_named_numeric_vector("543sp_Arth_Tpa_qval.txt")
Tge_qval_GL <- make_named_numeric_vector("543sp_Arth_Tge_qval.txt")
Tps_qval_GL <- make_named_numeric_vector("543sp_Arth_Tps_qval.txt")
Tdi_qval_GL <- make_named_numeric_vector("543sp_Arth_Tdi_qval.txt")


run_enrichment <- function(genelist, ref, sig_for_GO){
	
	### make rule for classing sig / non-sig 
	
	topDiffGenes <- function(allScore) {return(allScore < sig_for_GO)}
	
	#### make GOdata object
	#### setting node size as 10 so at least 10 genes must be annot per GO terms 
	#### do enrichment test
	
	GODATA_BP = new("topGOdata", ontology = "BP", allGenes = genelist, geneSel = topDiffGenes,  annot = annFUN.gene2GO, gene2GO = ref, nodeSize = 10)

	### get N GOs used

	GO_term_use_BP_list = GODATA_BP@graph@nodes
	N_GO_term_use_BP = length(GODATA_BP@graph@nodes)
    resultFisher <- runTest(GODATA_BP, algorithm = "weight01", statistic = "fisher")

	### combined tables
	allRes1_BP <- GenTable(GODATA_BP, Fisher_w01 = resultFisher, ranksOf = "Fisher_w01", topNodes = length(GODATA_BP@graph@nodes), numChar = 200)

	sig_fisher_BP_GO     = subset(allRes1_BP, allRes1_BP$Fisher_w01 < sig_for_GO)$GO.ID
		
	## return everything!
	out_list = list("N_GO_term_use_BP" = N_GO_term_use_BP, 
	                "GO_term_use_BP_list" = GO_term_use_BP_list, 
	                "allRes1_BP" = allRes1_BP, 
	                "sig_fisher_BP_GO" = sig_fisher_BP_GO,
	                "GODATA_BP" = GODATA_BP) 
	return(out_list)

}


#### run the enrichment stuff (0.05)

# # Tbi_Arth_enrich  <- run_enrichment(Tbi_qval_GL, geneID2GO_Tbi_Arth, 0.05)
# Tte_Arth_enrich  <- run_enrichment(Tte_qval_GL, geneID2GO_Tte_Arth, 0.05)
# Tce_Arth_enrich  <- run_enrichment(Tce_qval_GL, geneID2GO_Tce_Arth, 0.05)
# Tms_Arth_enrich  <- run_enrichment(Tms_qval_GL, geneID2GO_Tms_Arth, 0.05)
# Tcm_Arth_enrich  <- run_enrichment(Tcm_qval_GL, geneID2GO_Tcm_Arth, 0.05)
# Tsi_Arth_enrich  <- run_enrichment(Tsi_qval_GL, geneID2GO_Tsi_Arth, 0.05)
# Tpa_Arth_enrich  <- run_enrichment(Tpa_qval_GL, geneID2GO_Tpa_Arth, 0.05)
# Tge_Arth_enrich  <- run_enrichment(Tge_qval_GL, geneID2GO_Tge_Arth, 0.05)
# Tps_Arth_enrich  <- run_enrichment(Tps_qval_GL, geneID2GO_Tps_Arth, 0.05)
# Tdi_Arth_enrich  <- run_enrichment(Tdi_qval_GL, geneID2GO_Tdi_Arth, 0.05)

#### run the enrichment stuff (0.05)

Tbi_Droso_enrich  <- run_enrichment(Tbi_qval_GL, geneID2GO_Tbi_Droso, 0.05)
Tte_Droso_enrich  <- run_enrichment(Tte_qval_GL, geneID2GO_Tte_Droso, 0.05)
Tce_Droso_enrich  <- run_enrichment(Tce_qval_GL, geneID2GO_Tce_Droso, 0.05)
Tms_Droso_enrich  <- run_enrichment(Tms_qval_GL, geneID2GO_Tms_Droso, 0.05)
Tcm_Droso_enrich  <- run_enrichment(Tcm_qval_GL, geneID2GO_Tcm_Droso, 0.05)
Tsi_Droso_enrich  <- run_enrichment(Tsi_qval_GL, geneID2GO_Tsi_Droso, 0.05)
Tpa_Droso_enrich  <- run_enrichment(Tpa_qval_GL, geneID2GO_Tpa_Droso, 0.05)
Tge_Droso_enrich  <- run_enrichment(Tge_qval_GL, geneID2GO_Tge_Droso, 0.05)
Tps_Droso_enrich  <- run_enrichment(Tps_qval_GL, geneID2GO_Tps_Droso, 0.05)
Tdi_Droso_enrich  <- run_enrichment(Tdi_qval_GL, geneID2GO_Tdi_Droso, 0.05)


#####################################################################################################################
### Tidy


Tbi_Droso_enrich_table <- Tbi_Droso_enrich$allRes1_BP
Tte_Droso_enrich_table <- Tte_Droso_enrich$allRes1_BP
Tce_Droso_enrich_table <- Tce_Droso_enrich$allRes1_BP
Tms_Droso_enrich_table <- Tms_Droso_enrich$allRes1_BP
Tcm_Droso_enrich_table <- Tcm_Droso_enrich$allRes1_BP
Tsi_Droso_enrich_table <- Tsi_Droso_enrich$allRes1_BP
Tpa_Droso_enrich_table <- Tpa_Droso_enrich$allRes1_BP
Tge_Droso_enrich_table <- Tge_Droso_enrich$allRes1_BP
Tps_Droso_enrich_table <- Tps_Droso_enrich$allRes1_BP
Tdi_Droso_enrich_table <- Tdi_Droso_enrich$allRes1_BP

### add sp

Tbi_Droso_enrich_table$sp <- rep("Tbi", length(Tbi_Droso_enrich_table[,1]))
Tte_Droso_enrich_table$sp <- rep("Tte", length(Tte_Droso_enrich_table[,1]))
Tce_Droso_enrich_table$sp <- rep("Tce", length(Tce_Droso_enrich_table[,1]))
Tms_Droso_enrich_table$sp <- rep("Tms", length(Tms_Droso_enrich_table[,1]))
Tcm_Droso_enrich_table$sp <- rep("Tcm", length(Tcm_Droso_enrich_table[,1]))
Tsi_Droso_enrich_table$sp <- rep("Tsi", length(Tsi_Droso_enrich_table[,1]))
Tpa_Droso_enrich_table$sp <- rep("Tpa", length(Tpa_Droso_enrich_table[,1]))
Tge_Droso_enrich_table$sp <- rep("Tge", length(Tge_Droso_enrich_table[,1]))
Tps_Droso_enrich_table$sp <- rep("Tps", length(Tps_Droso_enrich_table[,1]))
Tdi_Droso_enrich_table$sp <- rep("Tdi", length(Tdi_Droso_enrich_table[,1]))

### join

All_Droso_enrich_table <- rbind(
Tbi_Droso_enrich_table,
Tte_Droso_enrich_table,
Tce_Droso_enrich_table,
Tms_Droso_enrich_table,
Tcm_Droso_enrich_table,
Tsi_Droso_enrich_table,
Tpa_Droso_enrich_table,
Tge_Droso_enrich_table,
Tps_Droso_enrich_table,
Tdi_Droso_enrich_table
)

### 

head(All_Droso_enrich_table)

All_Droso_enrich_table_sig <- subset(All_Droso_enrich_table, All_Droso_enrich_table$Fisher_w01 < 0.05)


All_Droso_enrich_table_sig %>% count(sp)

   # sp  n
# 1 Tbi  1
# 2 Tcm  4
# 3 Tdi  9
# 4 Tge  2
# 5 Tms  7
# 6 Tpa  2
# 7 Tsi  9
# 8 Tte 19


######################################################################
### N GO terms annot in sex and asex genes with sig +ve sel


get_sig_gene_vect <- function(list_file_name){
	df     <- read.table(list_file_name)
	df_sig <- subset(df, df$V2 < 0.05)
	return(as.character(df_sig$V1))
}

Tbi_sig_genes <- get_sig_gene_vect("543sp_Droso_Tbi_qval.txt")
Tte_sig_genes <- get_sig_gene_vect("543sp_Droso_Tte_qval.txt")
Tce_sig_genes <- get_sig_gene_vect("543sp_Droso_Tce_qval.txt")
Tms_sig_genes <- get_sig_gene_vect("543sp_Droso_Tms_qval.txt")
Tcm_sig_genes <- get_sig_gene_vect("543sp_Droso_Tcm_qval.txt")
Tsi_sig_genes <- get_sig_gene_vect("543sp_Droso_Tsi_qval.txt")
Tpa_sig_genes <- get_sig_gene_vect("543sp_Droso_Tpa_qval.txt")
Tge_sig_genes <- get_sig_gene_vect("543sp_Droso_Tge_qval.txt")
Tps_sig_genes <- get_sig_gene_vect("543sp_Droso_Tps_qval.txt")
Tdi_sig_genes <- get_sig_gene_vect("543sp_Droso_Tdi_qval.txt")



## get all BPs (i.e. not fitered by node size)
run_enrichment_for_all_BP <- function(genelist, ref, sig_for_GO){
	
	### make rule for classing sig / non-sig - note this rule is not used for the GSEA
	
	topDiffGenes <- function(allScore) {return(allScore < sig_for_GO)}
	# topDiffGenes <- function(allScore) {return(allScore < 1)} ## as a check - setting to one gives the same pvalues for the GSEA
	
	#### make GOdata object
	#### setting node size as 1 to get all BP gos 
	
	GODATA_BP = new("topGOdata", ontology = "BP", allGenes = genelist, geneSel = topDiffGenes,  annot = annFUN.gene2GO, gene2GO = ref, nodeSize = 1)
	### get N GOs used

	GO_term_use_BP_list = GODATA_BP@graph@nodes
	## return everything!
	out_list = list("GO_term_use_BP_list" = GO_term_use_BP_list, "GODATA_BP" = GODATA_BP) 
	return(out_list)

}


Tbi_Droso_all_BP  <- run_enrichment_for_all_BP(Tbi_qval_GL, geneID2GO_Tbi_Droso, 0.05)
Tte_Droso_all_BP  <- run_enrichment_for_all_BP(Tte_qval_GL, geneID2GO_Tte_Droso, 0.05)
Tce_Droso_all_BP  <- run_enrichment_for_all_BP(Tce_qval_GL, geneID2GO_Tce_Droso, 0.05)
Tms_Droso_all_BP  <- run_enrichment_for_all_BP(Tms_qval_GL, geneID2GO_Tms_Droso, 0.05)
Tcm_Droso_all_BP  <- run_enrichment_for_all_BP(Tcm_qval_GL, geneID2GO_Tcm_Droso, 0.05)
Tsi_Droso_all_BP  <- run_enrichment_for_all_BP(Tsi_qval_GL, geneID2GO_Tsi_Droso, 0.05)
Tpa_Droso_all_BP  <- run_enrichment_for_all_BP(Tpa_qval_GL, geneID2GO_Tpa_Droso, 0.05)
Tge_Droso_all_BP  <- run_enrichment_for_all_BP(Tge_qval_GL, geneID2GO_Tge_Droso, 0.05)
Tps_Droso_all_BP  <- run_enrichment_for_all_BP(Tps_qval_GL, geneID2GO_Tps_Droso, 0.05)
Tdi_Droso_all_BP  <- run_enrichment_for_all_BP(Tdi_qval_GL, geneID2GO_Tdi_Droso, 0.05)


get_N_GOs_wBPfilt <- function(want_vec_name, geneID2GO_file, BP_GOs_name, ALL_BP_GOs_name, sp, rep_m){
	out_df = c()
	for(i in seq(1:length(want_vec_name))){
		#print(want_vec_name[i])
		a1 <- eval(parse(text=paste(geneID2GO_file,'$',want_vec_name[i],sep='')))
		a1_filt <- a1[a1 %in% BP_GOs_name]
		a2_filt <- a1[a1 %in% ALL_BP_GOs_name]		
		#print(a1)
		#print(a1_filt)		
		#print("\n")
		out_df <- rbind(out_df, c(want_vec_name[i], length(a1), length(a2_filt), length(a1_filt)))
		colnames(out_df) <- c("gene_name", "N_All_GOs", "N_AllBP_GOs", "N_usedBP_GOs")
		} 
	out_df <- as.data.frame(out_df)
	out_df$N_All_GOs    <- as.numeric(as.character(out_df$N_All_GOs))
	out_df$N_AllBP_GOs  <- as.numeric(as.character(out_df$N_AllBP_GOs))
	out_df$N_usedBP_GOs <- as.numeric(as.character(out_df$N_usedBP_GOs))
	out_df$sp           <- rep(sp, length(out_df[,1]))	
	out_df$rep_mode     <- rep(rep_m, length(out_df[,1]))			
	return(out_df)
}


### get N GOs

Tbi_Droso_NGOs <- get_N_GOs_wBPfilt(Tbi_sig_genes, "geneID2GO_Tbi_Droso", Tbi_Droso_enrich$GO_term_use_BP_list, Tbi_Droso_all_BP$GO_term_use_BP_list, "Tbi", "sex")
Tte_Droso_NGOs <- get_N_GOs_wBPfilt(Tte_sig_genes, "geneID2GO_Tte_Droso", Tte_Droso_enrich$GO_term_use_BP_list, Tte_Droso_all_BP$GO_term_use_BP_list, "Tte", "asex")
Tce_Droso_NGOs <- get_N_GOs_wBPfilt(Tce_sig_genes, "geneID2GO_Tce_Droso", Tce_Droso_enrich$GO_term_use_BP_list, Tce_Droso_all_BP$GO_term_use_BP_list, "Tce", "sex")
Tms_Droso_NGOs <- get_N_GOs_wBPfilt(Tms_sig_genes, "geneID2GO_Tms_Droso", Tms_Droso_enrich$GO_term_use_BP_list, Tms_Droso_all_BP$GO_term_use_BP_list, "Tms", "asex")
Tcm_Droso_NGOs <- get_N_GOs_wBPfilt(Tcm_sig_genes, "geneID2GO_Tcm_Droso", Tcm_Droso_enrich$GO_term_use_BP_list, Tcm_Droso_all_BP$GO_term_use_BP_list, "Tcm", "sex")
Tsi_Droso_NGOs <- get_N_GOs_wBPfilt(Tsi_sig_genes, "geneID2GO_Tsi_Droso", Tsi_Droso_enrich$GO_term_use_BP_list, Tsi_Droso_all_BP$GO_term_use_BP_list, "Tsi", "asex")
Tpa_Droso_NGOs <- get_N_GOs_wBPfilt(Tpa_sig_genes, "geneID2GO_Tpa_Droso", Tpa_Droso_enrich$GO_term_use_BP_list, Tpa_Droso_all_BP$GO_term_use_BP_list, "Tpa", "sex")
Tge_Droso_NGOs <- get_N_GOs_wBPfilt(Tge_sig_genes, "geneID2GO_Tge_Droso", Tge_Droso_enrich$GO_term_use_BP_list, Tge_Droso_all_BP$GO_term_use_BP_list, "Tge", "asex")
Tps_Droso_NGOs <- get_N_GOs_wBPfilt(Tps_sig_genes, "geneID2GO_Tps_Droso", Tps_Droso_enrich$GO_term_use_BP_list, Tps_Droso_all_BP$GO_term_use_BP_list, "Tps", "sex")
Tdi_Droso_NGOs <- get_N_GOs_wBPfilt(Tdi_sig_genes, "geneID2GO_Tdi_Droso", Tdi_Droso_enrich$GO_term_use_BP_list, Tdi_Droso_all_BP$GO_term_use_BP_list, "Tdi", "asex")


Allsp_Droso_NGOs <- as.data.frame(rbind(
Tbi_Droso_NGOs,
Tte_Droso_NGOs,
Tce_Droso_NGOs,
Tms_Droso_NGOs,
Tcm_Droso_NGOs,
Tsi_Droso_NGOs,
Tpa_Droso_NGOs,
Tge_Droso_NGOs,
Tps_Droso_NGOs,
Tdi_Droso_NGOs
))


Allsp_Droso_NGOs$N_All_GOs_bi <- ifelse(Allsp_Droso_NGOs$N_All_GOs > 0, 1,0)
Allsp_Droso_NGOs$N_AllBP_GOs_bi <- ifelse(Allsp_Droso_NGOs$N_AllBP_GOs > 0, 1,0)
Allsp_Droso_NGOs$N_usedBP_GOs_bi <- ifelse(Allsp_Droso_NGOs$N_usedBP_GOs > 0, 1,0)


Allsp_Droso_NGOs$sp_pair <-
ifelse(Allsp_Droso_NGOs$sp == "Tbi", "Tbi-Tte",
ifelse(Allsp_Droso_NGOs$sp == "Tce", "Tce-Tms",
ifelse(Allsp_Droso_NGOs$sp == "Tcm", "Tcm-Tsi",
ifelse(Allsp_Droso_NGOs$sp == "Tpa", "Tpa-Tge",
ifelse(Allsp_Droso_NGOs$sp == "Tps", "Tps-Tdi",
ifelse(Allsp_Droso_NGOs$sp == "Tte", "Tbi-Tte", 
ifelse(Allsp_Droso_NGOs$sp == "Tms", "Tce-Tms",
ifelse(Allsp_Droso_NGOs$sp == "Tsi", "Tcm-Tsi",
ifelse(Allsp_Droso_NGOs$sp == "Tge", "Tpa-Tge",
ifelse(Allsp_Droso_NGOs$sp == "Tdi", "Tps-Tdi", 
NA))))))))))

Allsp_Droso_NGOs$sp_pair_ord = ordered(Allsp_Droso_NGOs$sp_pair, levels = c("Tbi-Tte", "Tcm-Tsi", "Tce-Tms", "Tps-Tdi", "Tpa-Tge"))



head(Allsp_Droso_NGOs)

### at least one GO term annot 
Allsp_Droso_NGOs_more_than_0_used <- subset(Allsp_Droso_NGOs, Allsp_Droso_NGOs$N_usedBP_GOs > 0)


### plot

P1_N_usedBP_GOs_prop_GOs <- ggplot(Allsp_Droso_NGOs)  + 
	geom_bar(aes(sp_pair_ord, N_usedBP_GOs_bi, fill = as.factor(rep_mode)), position = "dodge", stat = "summary", fun = "mean") +
	theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
	scale_fill_manual(values = c("#92C5DECD", "#D6604DED")) + 
	xlab ("Species pair") + 
	ylab ("Prop genes with GO terms annotated")  

P1_N_usedBP_GOs_mean_GOs_more_than_0_BPused_GO <-  ggplot(Allsp_Droso_NGOs_more_than_0_used )  + 
	geom_bar(aes(sp_pair_ord, N_usedBP_GOs, fill = as.factor(rep_mode)), position = "dodge", stat = "summary", fun = "mean") +
	theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
	scale_fill_manual(values = c("#92C5DECD", "#D6604DED")) + 
	xlab ("Species pair") + 
	ylab ("Mean GO terms annotated per gene")  


plot_grid(P1_N_usedBP_GOs_prop_GOs, P1_N_usedBP_GOs_mean_GOs_more_than_0_BPused_GO)























