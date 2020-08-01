library('RColorBrewer')
library('AsexStats')

pdf('figures/SNP_pairwise_similarities.pdf')

sp <- '5_Tge'

tab_filename <- paste0('data/SNP_calls/', sp ,'_reduced_variants.tsv')
variant_tab <- read.table(tab_filename, stringsAsFactors = F)
colnames(variant_tab) <- c('scf', 'pos', 'qual', paste0('g', 1:5), paste0('d', 1:5))

variant_tab$g0 <- '0/0'

genotypes <- c('g0', 'g1', 'g2', 'g3', 'g4', 'g5')

similarity_matrix <- matrix(NA, nrow = 6, ncol = 6)
dimnames(similarity_matrix) <- list(genotypes, genotypes)

shared_genotypes <- function(g_vec1, g_vec2){
    filter_out <- g_vec1 == './.' | g_vec2 == './.'
    sum(g_vec1[!filter_out] == g_vec2[!filter_out])
}

for ( ind_1 in genotypes ){
    for ( ind_2 in genotypes ){
        if ( is.na(similarity_matrix[ind_2, ind_1]) ){
            similarity_matrix[ind_1, ind_2] <- shared_genotypes(variant_tab[, ind_1], variant_tab[, ind_2])
        } else {
            similarity_matrix[ind_1, ind_2] <- similarity_matrix[ind_2, ind_1]
        }
    }
}

individuals <- paste0(substr(sp, 3, 5), '_0', 0:5)
dimnames(similarity_matrix) <- list(individuals, individuals)

heatmap(similarity_matrix, main = sp)