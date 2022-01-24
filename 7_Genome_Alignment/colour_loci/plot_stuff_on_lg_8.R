library(AsexStats)
library(RColorBrewer)

window = 1e5
gap_beween_chromosomes = 3

plot_SVs_on_LGs <- function(exclude_zero = T, lines = F, sex = 'asex', SNP_ylim, SV_ylim, xlim = NA){

    if ( is.na(xlim) ){
        xlim <- c(1, nrow(variant_density_table))
    }

    title <- timemas$labels[sp == timemas$codes]

    if ( sex == 'asex' ){

        pal <- brewer.pal(5, "YlGnBu")[c(3,5)]
        plot(NULL, xlim = xlim, ylim = c(0, 1), pch = 20,
             main = title,
             xaxt = "n", yaxt = "n", bty = 'n', xlab = '', ylab = '', cex.axis = 1.4, cex.main = 1.6)
    } else {
        pal <- brewer.pal(5, "YlOrRd")[c(3,5)]
        plot(NULL, xlim = xlim, ylim = c(0, 1), pch = 20,
             xaxt = "n", yaxt = "n", bty = 'n', xlab = '', ylab = '', cex.axis = 1.4, cex.main = 1.6)
    }

     # main = timemas$labels[timemas$codes == sp],
     # xlab = 'linage group [ Mbp ]', ylab = '# found SVs'
    xtick <- chromosomes$adjustments / window
    # axis(side = 1, at = xtick, labels = FALSE)
    if ( sp == '5_Tpa'){
        text(x = (chromosomes$adjustments[1:12] + chromosomes$adjustments[2:13]) / (2 * window),  par("usr")[3],
             labels = 1:12, pos = 1, xpd = TRUE, cex = 1, xpd=NA)
        text(x = (chromosomes$adjustments[11] + chromosomes$adjustments[12]) / (2 * window) + 50,  par("usr")[3],
             labels = 12, pos = 1, xpd = TRUE, cex = 1, xpd=NA)
    }

    row.names(ref_scf2lg_pos) <- ref_scf2lg_pos$scf
    rect((ref_scf2lg_pos['lg8_ord14_scaf702.1', 'to'] / window) - 41, -2, (ref_scf2lg_pos['lg8_ord15_scaf128', 'from'] / window) + 60, 1, col = 'lightgrey', border = F)
    # if ( sp == '2_Tcm'){
    rect((ref_scf2lg_pos['lg8_ord15_scaf128', 'from'] / window) + 50, -2, (ref_scf2lg_pos['lg8_ord15_scaf128', 'from'] / window) + 60, 1, col = 'khaki', border = F)
    # } else {
    rect((ref_scf2lg_pos['lg8_ord15_scaf128', 'from'] / window) + 51.6, -2, (ref_scf2lg_pos['lg8_ord15_scaf128', 'from'] / window) + 56.69, 1, col = 'cornsilk', border = F)
    # }

    # if ( sex == 'asex'){
        axis(side = 2, at = seq(0,1, by = 0.2), labels = seq(0, SNP_ylim, length = 6))
        mtext('SNPs', 2, line = 2.5)
    # } else {
        axis(side = 4, at = seq(0,1, by = 0.2), labels = seq(0, SV_ylim, length = 6))
        mtext('SVs', 4, line = 2.5)
    # }

    if ( exclude_zero ){
        variant_density_table[variant_density_table$SNPs == 0, 'SNPs'] <- NA
        variant_density_table[variant_density_table$SVs == 0, 'SVs'] <- NA
    }

    if ( lines ){
        lines(variant_density_table$SVs / SV_ylim, col = pal[2], lwd = 1.6)
        lines(variant_density_table$SNPs / SNP_ylim, col = pal[1], lwd = 1.6)
    } else {
        points(variant_density_table$SVs / SV_ylim, pch = 20, col = pal[2])
        points(variant_density_table$SNPs / SNP_ylim, pch = 20, col = pal[1])
    }
    # legend('topright', bty = 'n', c('SNPs', 'SVs'), pch = 20, col = pal[c(1,2)])
}


# interesting species: 2_Tcm; 3_Tms
# potentially also 4_Tbi, 4_Tte, 5_Tpa
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

    # keep only lg8 blocks
    mapping_table <- mapping_table[mapping_table$lg == 'lg8', ]
    nrow(mapping_table)

    mapping_table <- mapping_table[order(mapping_table$lg, mapping_table$lg_start), ]

    scaf128 <- mapping_table[mapping_table$ref == "lg8_ord15_scaf128",]

    scaf128[scaf128$block_r_start < 5.69e6 & scaf128$block_r_end > 5.16e6, ]
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
            if (lg == 'lg8'){
                lg_pos <- sv_frame[i, 'pos'] + relevant_subset[mapped_to, 'lg_start']
                row <- variant_density_table$lg == lg & variant_density_table$lg_from < lg_pos & variant_density_table$lg_to > lg_pos
                variant_density_table[row, 'SVs'] = variant_density_table[row, 'SVs'] + 1
            }
        }
        if (sum(mapped_to) > 1) {
            multimapped_filt = multimapped_filt + 1
        }
    }

    print(paste('Removing', multimapped_filt, 'variants due to multimapping'))

    ### write the table down

    write.table(variant_density_table, output_file, quote = F, sep = '\t')

    plot_SVs_on_LGs(F, F, 'sex', 6000, 25, c(450.43989, 550.43989))

    for(i in 1:nrow(ref_scf2lg_pos)){
        lines(c(ref_scf2lg_pos$from[i] / window, ref_scf2lg_pos$to[i] / window), c(0,0), lwd = 3, col = 'red')
        lines(c(ref_scf2lg_pos$from[i] / window, ref_scf2lg_pos$from[i] / window), c(-0.02,0.02), lwd = 3, col = 'red')
        lines(c(ref_scf2lg_pos$to[i] / window, ref_scf2lg_pos$to[i] / window), c(-0.02,0.02), lwd = 3, col = 'red')
    }

    label_pos <- ((ref_scf2lg_pos$from / window) + (ref_scf2lg_pos$to / window)) / 2
    labels <- sapply(strsplit(ref_scf2lg_pos$scf, '_'), function(x){ x[3] })
    # axis(side = 1, at = xtick, labels = FALSE)
    text(x = label_pos,  par("usr")[3], labels = labels, pos = 1, xpd = TRUE, cex = 1, xpd=NA, srt = 90)

}