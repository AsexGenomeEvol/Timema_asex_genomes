str2depth <- function(genotype_str, index = 4, which_genotype = 1){
    if ( is.na(genotype_str[index]) ){
        return(NA)
    } else {
        as.numeric(unlist(strsplit(genotype_str[index], ','))[which_genotype])
    }
}

plot_one <- function(SV_tab, type, main){
    SV_tab$genotype <- substr(SV_tab[,10], 1, 3)
    element <- ifelse(type == 'SR', 6, 5)
    SV_tab$cov_1 <- sapply(strsplit(SV_tab[,10], ':'), str2depth, element, 1)
    SV_tab$cov_2 <- sapply(strsplit(SV_tab[,10], ':'), str2depth, element, 2)
    SV_tab$cov <- SV_tab$cov_1 + SV_tab$cov_2

    ggplot(SV_tab, aes(cov, fill = genotype)) + geom_density(alpha = 0.2) + xlim(0, 100) + ggtitle(main)
}