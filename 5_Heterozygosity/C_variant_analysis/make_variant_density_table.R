library(AsexStats)

window = 1e6
gap_beween_chromosomes = 3

# sp = '3_Tms'
for(sp in timemas$codes){

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
    nrow(mapping_table)

    # filter out blocks mapped "unequally" - reference block / query block can not be < 0.5
    mapping_table <- mapping_table[abs(mapping_table$block_r_end - mapping_table$block_r_start) / mapping_table$block_size > 0.5, ]
    nrow(mapping_table)

    mapping_table <- mapping_table[order(mapping_table$lg, mapping_table$lg_start), ]

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

    variant_density_table <- do.call("rbind", lapply(1:12, get_lg_windows))

    ########################
    ### add mapping info ###
    ########################

    variant_density_table$uniq_mapped <- 0
    lg_mapping_table <- data.frame(lg = paste0('lg', 1:12), unique_mapped = NA, multimapped = NA)
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

    variant_density_table$uniq_mapped[is.na(variant_density_table$uniq_mapped)] <- 0
    variant_density_table$uniq_mapped <- round(variant_density_table$uniq_mapped, 3)

    ####################
    ### add SNP info ###
    ####################

    tab_filename <- paste0('data/SNP_calls/', sp, '_reduced_filtered_variants.tsv')

    variant_tab <- read.table(tab_filename, stringsAsFactors = F)
    colnames(variant_tab) <- c('scf', 'pos', 'qual', paste0('g', 1:5), paste0('d', 1:5), 'ref_scf', 'ref_pos', 'lg', 'lg_pos')

    if ( sp == '1_Tps'){
        variant_tab[,'g4'] <- NA
        variant_tab[,'g5'] <- NA
    }
    if ( sp == '4_Tte'){
        variant_tab[,'g1'] <- NA
    }
    if ( sp == '2_Tsi'){
        variant_tab[,'g1'] <- NA
        variant_tab[,'g2'] <- NA
        variant_tab[,'g3'] <- NA
    }

    print(paste('Filtering ', round(100 * mean(is.na(variant_tab$lg_pos)), 1), ' % of variants (unmapped)'))
    anchored_variant_tab <- variant_tab[!is.na(variant_tab$lg_pos), ]

    variant_density_table$SNPs <- 0

    for (i in 1:nrow(anchored_variant_tab)) {
        lg = anchored_variant_tab[i, 'lg']
        pos = anchored_variant_tab[i, 'lg_pos']
        row <- variant_density_table$lg == lg & variant_density_table$lg_from < pos & variant_density_table$lg_to > pos
        variant_density_table[row, 'SNPs'] <- variant_density_table[row, 'SNPs'] + 1
    }

    ###################
    ### add SV info ###
    ###################

    source('D_variant_analysis/load_SV_calls.R')

    SVs_filt_stringent_union_file <- paste0('data/manta_SV_calls/data/', sp, '/SVs_filt_stringent_union.vcf')

    one_sp_sv_calls <- load_SV_calls(SVs_filt_stringent_union_file)[[1]]
    info_vec <- sapply(one_sp_sv_calls, function(x){x[8]} )
    info_vec <- strsplit(info_vec, ';')

    sv_frame <- data.frame(
        scf = sapply(info_vec, get_entry, "CHR2"),
        pos = sapply(one_sp_sv_calls, function(x) { as.numeric(x[2]) } ),
        len = abs(as.numeric(sapply(info_vec, get_entry, "SVLEN"))),
        type = sapply(info_vec, get_entry, "SVTYPE")
    )

    presence_vec = sapply(one_sp_sv_calls, get_genotypes)
    sv_frame[, c('01', '02', '03', '04', '05')] <- matrix(presence_vec, ncol = 5, byrow = T)

    # I need to calculate number of anchored SVs per window
    variant_density_table$SVs <- 0
    multimapped_filt = 0

    for (i in 1:nrow(sv_frame)) {
        # find it in the mapping_table
        relevant_subset <- mapping_table[mapping_table$scf == sv_frame[i, 'scf'], ]
        sv_from <- sv_frame[i, 'pos']
        sv_to <- sv_frame[i, 'pos'] + sv_frame[i, 'len']
        block_froms <- apply(relevant_subset[, c('block_q_start', 'block_q_end')], 1, min)
        block_tos <- apply(relevant_subset[, c('block_q_start', 'block_q_end')], 1, max)

        mapped_to <- (block_froms < sv_from & sv_from < block_tos | block_froms < sv_to & sv_to < block_tos)

        if (sum(mapped_to) == 1){
            # unique mapping
            lg <- relevant_subset[mapped_to, 'lg']
            lg_pos <- sv_frame[i, 'pos'] + relevant_subset[mapped_to, 'lg_start']
            row <- variant_density_table$lg == lg & variant_density_table$lg_from < lg_pos & variant_density_table$lg_to > lg_pos
            variant_density_table[row, 'SVs'] = variant_density_table[row, 'SVs'] + 1
        }
        if (sum(mapped_to) > 1) {
            multimapped_filt = multimapped_filt + 1
        }
    }

    print(paste('Removing', multimapped_filt, 'variants due to multimapping'))

    ### write the table down

    write.table(variant_density_table, output_file, quote = F, sep = '\t')

}