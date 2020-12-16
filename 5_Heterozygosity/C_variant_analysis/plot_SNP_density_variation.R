library(AsexStats)
library(DescTools)

window = 1e6

variant_density_table_files <- paste0("tables/", timemas$codes, "_variants_on_chromosomes_w", window, ".tsv")
variant_density_tables <- lapply(variant_density_table_files, read.table, header = T, sep = '\t', stringsAsFactors = F)

get_mad <- function(one_tab, method = 'all', ci = NA){
    snp_densities <- one_tab$SNPs / (one_tab$uniq_mapped * 1e6)
    snp_densities <- snp_densities[!is.na(snp_densities)]
    snp_densities <- snp_densities[is.finite(snp_densities)]

    if ( method == 'all'){
        if ( is.na(ci) ){
            return(mad(snp_densities))
        }
        return(MADCI(snp_densities))
    }
}

variances <- data.frame(sp = timemas$codes)
variances[, c('mad', 'lwr.ci', 'upr.ci')] <-  t(sapply(variant_density_tables, get_mad, 'all', T))

ymax <- max(variances[, c('upr.ci')])
locations <- barplot(variances$mad, col = c(asex_blue, sex_red), ylab = 'Median Absolute Deviation', ylim = c(0, ymax))
text(locations,
     par("usr")[3] - (0.03 * ymax), pos = 1, srt = 25,
     xpd = TRUE, labels = timemas$labels, cex = 0.8)

width = c(-0.1, 0.1)

for ( i in 1:10 ){
    lines(c(locations[i], locations[i]), c(variances[i, 'lwr.ci'], variances[i, 'upr.ci'] ))
    lines(c(locations[i], locations[i]) + width, c(variances[i, 'lwr.ci'], variances[i, 'lwr.ci'] ))
    lines(c(locations[i], locations[i]) + width, c(variances[i, 'upr.ci'], variances[i, 'upr.ci'] ))
}

# variances$pair <- substr(variances$sp, 1, 1)
# variances$repr_mode <- c('asex', 'sex')
#
# repr_mode_test <- glm(mad ~ pair + repr_mode, data = variances)

mad_asex <- variances$mad[seq(1,10, by = 2)]
mad_sex <- variances$mad[seq(2,10, by = 2)]

repr_mode_test <- t.test(mad_sex, mad_asex, paired = TRUE, alternative = "greater")
summary(repr_mode_test)