ref_ind_file <- c('variant_calls/Tms_00/Tms_00_union.summary',
                  'variant_calls/Tms_00/Tms_00_merged.summary')

asm_errrors_file <- c('variant_calls/delly/merged_all.summary',
                      'variant_calls/lumpy/merged_all.summary',
                      'variant_calls/manta/merged_all.summary')

asm_errror_union_file <- ('asm_errors_union.summary')

read.call.summary <- function(file){
    read.table(file, sep = '\t', header = T)
}

asm_err_union <- read.call.summary(asm_errror_union_file)
ref_calls  <- lapply(ref_ind_file, read.call.summary)
asm_err_calls <- lapply(asm_errrors_file, read.call.summary)


sum(as.matrix(ref_calls[[1]][, c(2,3,4,5,6,7)]))
sum(as.matrix(ref_calls[[2]][, c(2,3,4,5,6,7)]))
sum(as.matrix(asm_err_union[, c(2,3,4,5,6,7)]))

delly_merged <- read.table("variant_calls/delly/merged_merged.SUPP_VEC", header = F, colClasses = 'character', stringsAsFactors = F)$V1
manta_merged <- read.table("variant_calls/manta/merged_merged.SUPP_VEC", header = F, colClasses = 'character', stringsAsFactors = F)$V1
lumpy_merged <- read.table("variant_calls/lumpy/merged_merged.SUPP_VEC", header = F, colClasses = 'character', stringsAsFactors = F)$V1

table(delly_merged)
barplot(table(delly_merged))

table(manta_merged)
barplot(table(manta_merged))

table(lumpy_merged)
barplot(table(lumpy_merged))


delly_matrix <- do.call(rbind, lapply(strsplit(delly_merged, split = ''), as.numeric))

####
manta_merged <- read.table("manta_merged.SUPP_VEC", header = F, colClasses = 'character', stringsAsFactors = F)$V1
lumpy_merged <- read.table("lumpy_merged.SUPP_VEC", header = F, colClasses = 'character', stringsAsFactors = F)$V1

table(manta_merged)
barplot(table(manta_merged))

table(lumpy_merged)
barplot(table(lumpy_merged))

##### Tce
menta_calls <- readLines("data/3_Tce/variant_calls/manta_merged.vcf")
menta_calls <- menta_calls[!grepl("^##", menta_calls)]
header <- unlist(strsplit(menta_calls[1], '\t'))

menta_calls <- strsplit(menta_calls[-1], '\t')
to_keep <- which(sapply(menta_calls, function(x) {grepl("SUPP_VEC=111111", x[8])} ))
menta_calls_in_all <- menta_calls[to_keep]

fractions_of_het <- sapply(menta_calls_in_all, function(x) { sum(grepl("0/1", x[10:15])) } )

library(AsexStats)
png('~/Desktop/Tce_number_of_heterozygots_almong_variants_shared_by_all.png')
    barplot(table(fractions_of_het), main = "Tce number of heterozygots almong variants shared by all")
dev.off()

alt_alleles <- function(x){
    sum(grepl("0/1", x[10:15]) * 1) +
    sum(grepl("1/1", x[10:15]) * 2)
}

site_freq <- sapply(menta_calls, alt_alleles)

png('~/Desktop/Tce_site_freq_spectra.png')
    barplot(table(site_freq), main = 'Tce site freq spectra (all SVs, no size filt, only SVs in at least 2 ind retained)', col = sex_red)
dev.off()


gen_heterozygots <- function(x){
    sum(grepl("0/1", x[10:15]))
}

pop_heterozygosity <- sapply(menta_calls, gen_heterozygots)

library(MASS)
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
cells <- table(paste(pop_heterozygosity, site_freq))
cell_coordinates <- lapply(strsplit(names(cells), ' '), as.numeric)

k <- kde2d(pop_heterozygosity, site_freq, n=100)

png('~/Desktop/Tms_alle_freq_vs_heterozyg.png')
    image(k, col=rf(100)[1:80],
          main = 'Tms heterozygosites vs allele freq', xlab = 'number of heterozygous individuals', ylab = 'allele frequency',
          xlim = c(-0.1, 6.1), ylim = c(1.8, 12.2))
    text(sapply(cell_coordinates, function(x){ if(x[1] == 0){ 0.05 } else { x[1]}} ), sapply(cell_coordinates, function(x){ x[2] } ), cells)
dev.off()

# like this I could substitute patterns in names
# test = c("asd", 'vfd')
# sapply(c("1", "2"), function(x) { paste(test, collapse = x) } )
