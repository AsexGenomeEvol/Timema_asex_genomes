library(AsexStats)
library(RColorBrewer)

window = 2e5
gap_beween_chromosomes = 3

sp <- '2_Tcm'
output_file <- paste0("tables/", sp, "_variants_on_chromosomes_w", window, ".tsv")

if ( file.exists(output_file) ){
    print(paste('The file:', output_file, 'exists already'))
    next
} else {
    print(paste('Generating ', output_file, 'file'))
}

# using the mapping ...
mapping_file <- paste0('data/b3v08_anchoring_to_LGs/', sp, '_scf_block_alignment.tsv')
mapping_table <- read.table(mapping_file, header = T)
mapping_table <- mapping_table[abs(mapping_table$block_r_end - mapping_table$block_r_start) / mapping_table$block_size > 0.5, ]
mapping_table <- mapping_table[mapping_table$lg == 'lg8', ]
mapping_table <- mapping_table[order(mapping_table$lg, mapping_table$lg_start), ]

# scaf128 <- mapping_table[mapping_table$ref == "lg8_ord15_scaf128",]
# scaf128[scaf128$block_r_start < 5.69e6 & scaf128$block_r_end > 5.16e6, ]
### generate blank variant_density_table

get_lg_windows <- function(i) {
    all_windows <- seq(0, chromosomes[i,'rounded_len'], by = window)
    adj <- chromosomes[i, 'adjustments']
    lg_from <- all_windows[1:(length(all_windows) - 1)]
    lg_to <- all_windows[2:length(all_windows)]
    data.frame(lg = chromosomes[i, 'chr'],
               lg_from = lg_from,
               lg_to = lg_to,
               genome_from = lg_from + adj,
               genome_to = lg_from + adj)
}

reference <- read.table('data/external_ref/sex_lg_assigment_scores_1.4a.tsv')
colnames(reference) <- c('scf_o', 'scf', 'score', 'cov', 'len', 'asignment')
reference$chromosome <- sapply(strsplit(reference$scf, "_"), function(x) { x[1] } )

source('D_variant_analysis/load_chromosomes.R')
chromosomes$adjustments <- 0

variant_density_table <- get_lg_windows(8)

########################
### add mapping info ###
########################

variant_density_table$uniq_mapped <- 0
lg_mapping_table <- data.frame(lg = 'lg8', unique_mapped = NA, multimapped = NA)
for (lg in lg_mapping_table$lg){
    one_lg <- mapping_table[mapping_table$lg == lg, ]
    one_lg <- one_lg[order(one_lg$lg_start), ]
    lg_nts <- rep(0, max(one_lg$lg_end))
    for (i in 1:nrow(one_lg)){
        lg_nts[one_lg$lg_start[i]:one_lg$lg_end[i]] <- lg_nts[one_lg$lg_start[i]:one_lg$lg_end[i]] + 1
    }
    lg_mapping_table[lg_mapping_table$lg == lg, 'unique_mapped'] <- mean(lg_nts == 1)
    lg_mapping_table[lg_mapping_table$lg == lg, 'multimapped'] <- mean(lg_nts > 1)
    for ( lg_row in which(variant_density_table$lg == lg) ){
        variant_density_table[lg_row, 'uniq_mapped'] <- mean(lg_nts[(variant_density_table[lg_row, 'lg_from'] + 1):variant_density_table[lg_row, 'lg_to']] == 1, rm.na = T)
    }

}

variant_density_table <- variant_density_table[!is.na(variant_density_table$uniq_mapped),]
variant_density_table$uniq_mapped <- round(variant_density_table$uniq_mapped, 3)

####################
### add SNP info ###
####################

tab_filename <- paste0('data/SNP_calls/', sp, '_reduced_filtered_variants.tsv')

variant_tab <- read.table(tab_filename, stringsAsFactors = F)
colnames(variant_tab) <- c('scf', 'pos', 'qual', paste0('g', 1:5), paste0('d', 1:5), 'ref_scf', 'ref_pos', 'lg', 'lg_pos')

print(paste('Filtering ', round(100 * mean(is.na(variant_tab$lg_pos)), 1), ' % of variants (unmapped)'))
anchored_variant_tab <- variant_tab[!is.na(variant_tab$lg_pos), ]

print(paste('Filtering ', round(100 * mean(anchored_variant_tab$lg != 'lg8'), 1), ' % of mapped variants (not mapped to lg8)'))
anchored_variant_tab <- anchored_variant_tab[anchored_variant_tab$lg == 'lg8', ]

variant_density_table$SNPs <- 0

for (i in 1:nrow(anchored_variant_tab)) {
    lg = anchored_variant_tab[i, 'lg']
    pos = anchored_variant_tab[i, 'lg_pos']
    row <- variant_density_table$lg == lg & variant_density_table$lg_from < pos & variant_density_table$lg_to > pos
    variant_density_table[row, 'SNPs'] <- variant_density_table[row, 'SNPs'] + 1
}

# getting reference scaffolds -> lg info
ref_scf2lg_pos = data.frame(scf = unique(anchored_variant_tab$ref_scf))
ref_scf2lg_pos$from <- NA
ref_scf2lg_pos$to <- NA

for (scf in ref_scf2lg_pos$scf){
    one_scf <- anchored_variant_tab[anchored_variant_tab$ref_scf == scf, ]
    ref_scf2lg_pos[scf == ref_scf2lg_pos$scf, 'from'] <- min(one_scf$lg_pos)
    ref_scf2lg_pos[scf == ref_scf2lg_pos$scf, 'to']   <- max(one_scf$lg_pos)
}

ref_scf2lg_pos <- ref_scf2lg_pos[order(ref_scf2lg_pos$from), ]

### write the table down

head(variant_density_table)

library('changepoint')

plot(variant_density_table$SNPs)

plot(cpt.mean(variant_density_table$SNPs, method = "BinSeg", minseglen = 80))
cpt.mean(variant_density_table$SNPs, method = "BinSeg", minseglen = 80)

