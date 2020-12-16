library('RColorBrewer')
library('AsexStats')

pdf('figures/SNP_pairwise_similarities.pdf')

asexuals <- timemas$codes[seq(1, 10, by=2)]

for (sp in timemas$codes){
    tab_filename <- paste0('data/SNP_calls/', sp, '_reduced_filtered_variants.tsv')
    variant_tab <- read.table(tab_filename, stringsAsFactors = F)
    colnames(variant_tab) <- c('scf', 'pos', 'qual', paste0('g', 1:5), paste0('d', 1:5))

    variant_tab$g0 <- '0/0'

    if ( sp == '1_Tps'){
        genotypes <- c('g0', 'g1', 'g2', 'g3')
    } else if ( sp == '4_Tte'){
        genotypes <- c('g0', 'g2', 'g3', 'g4', 'g5')
    } else if ( sp == '2_Tsi'){
        genotypes <- c('g0', 'g4', 'g5')
    } else {
        genotypes <- c('g0', 'g1', 'g2', 'g3', 'g4', 'g5')
    }

    different_genotypes <- function(g_vec1, g_vec2){
        filter_out <- g_vec1 == './.' | g_vec2 == './.'
        sum(g_vec1[!filter_out] != g_vec2[!filter_out])
    }

    variant_tab$alt_alleles <- apply(variant_tab[, genotypes], 1, function(x){ sum(x == '1/1') } )
    common_variant_tab <- variant_tab[variant_tab$alt_alleles %in% c(2,3,4), ]

    dist_matrix <- matrix(NA, nrow = length(genotypes), ncol = length(genotypes))
    dimnames(dist_matrix) <- list(genotypes, genotypes)
    for ( ind_1 in genotypes ){
        for ( ind_2 in genotypes ){
            if ( is.na(dist_matrix[ind_2, ind_1]) ){
                dist_matrix[ind_1, ind_2] <- different_genotypes(common_variant_tab[, ind_1], common_variant_tab[, ind_2])
            } else {
                dist_matrix[ind_1, ind_2] <- dist_matrix[ind_2, ind_1]
            }
        }
    }

    # dentrogram
    dentrogram <- hclust(as.dist(dist_matrix), method = "single")
    plot(dentrogram)

    # individuals <- paste0(substr(sp, 3, 5), '_0', substr(genotypes, 2, 2))
    # dimnames(dist_matrix) <- list(individuals, individuals)
    # colors_in_gradient <- 100
    # resolution <- seq(0, 1250000, len = colors_in_gradient)
    # pal <- colorRampPalette(c('white', asex_blue))(colors_in_gradient - 1)
    # genotype_matrix <- as.matrix(variant_tab[,genotypes])
    # transform for plotting
    # dist_matrix <- t(dist_matrix)[rev(1:nrow(dist_matrix)),]
    # image(, col = pal, main = sp, axes = F)

    # heatmap(mat, Colv=NA, col=greenred(10), scale=”none”)

}

get_derived_variant <- function(x){
    if (sum(x == '1/1') >= 3 ){
        paste(genotypes[which(x == '0/0')], collapse = "_")
    } else {
        paste(genotypes[which(x == '1/1')], collapse = "_")
    }
}

shared_derived_alleles <- apply(common_variant_tab[, genotypes], 1, get_derived_variant)

sort(table(shared_derived_alleles), decreasing = T)
# g0_g3_g4    g0_g5    g4_g5 g0_g2_g5    g0_g2    g3_g4    g1_g2    g1_g3
#   347376    64107    42091    36850    34234    31026    27675    23862
#    g2_g5 g0_g1_g2 g0_g2_g4    g3_g5    g2_g4    g2_g3 g0_g1_g5 g0_g4_g5
#    18821    15534    14763    13982    13374    12549    12254    11998
#    g0_g4    g1_g5    g1_g4 g0_g2_g3 g0_g3_g5 g0_g1_g3    g0_g3 g0_g1_g4
#    11850    11850     9557     8701     8334     6644     6387     5820
#    g0_g1
#     5709

genotypes <- c('g0', 'g1', 'g2', 'g3', 'g4')
shared_derived_alleles <- apply(common_variant_tab[, genotypes], 1, get_derived_variant)

sort(table(shared_derived_alleles), decreasing = T)
#  g1_g2  g0_g2  g3_g4  g1_g3  g0_g4  g2_g4  g2_g3  g1_g4
# 375319  70962  46642  38707  23700  20054  18424  18382
#  g0_g1  g0_g3
#  17927  14565


# paste(genotypes[which(common_variant_tab[1, genotypes] == '1/1')], collapse = "_")



dev.off()