library(AsexStats)

# load species genome info
scf_len_files <- paste0("data/", timemas$codes, "/reference/", timemas$codes, "_b3v08_scf.lengths")
scf_lengths <- lapply(scf_len_files, read.table, header = F, row.names = 1, col.names = c('scf', 'len'))
scf_to_analyze <- lapply(scf_lengths, function(x){rownames(x)[x[,1] > 100000]})

########
# load SV tables
source('F_structural_variation_analysis/load_SV_calls.R')
source('F_structural_variation_analysis/filter_SV_calls.R')

SV_call_manta_files <- paste0("data/", timemas$codes, "/variant_calls/", timemas$codes, "_survivor_manta_calls_union.vcf")
SV_calls <- load_SV_calls(SV_call_manta_files)

sv_genotype_tables <- lapply(SV_calls, get_sv_table, genotypes = T, T)

######
# args <- commandArgs(trailingOnly=TRUE)
sp_i <- 8

for (sp_i in 1:10){
    sp <- timemas$codes[sp_i]
    # parameters
    treshold_for_chaining <- 10000
    cov_treshold <- 0.6

    source('../timema_assembly/scripts/prepare_reference.R')
    # ??

    coord1_files <- paste0('data/', timemas$codes, '/genome_aln/to_ref/',
                           timemas$codes, '_b3v08_to_ref.1coords')
    output_file <- paste0("data/", sp, "/variant_calls/",
                          "3_Tms_survivor_manta_calls_union_mapped_on_LGs.tsv")

    mummer <- load_mummer_coords(file(coord1_files[sp_i]))
    # reducing the size of the queried dataset
    mummer <- mummer[mummer$q_name %in% scf_to_analyze[[sp_i]],]

    SV_table <- sv_genotype_tables[[sp_i]]
    # filtering out SVs found only in a single individual
    SV_table <- SV_table[rowSums(SV_table[,c('00','01','02','03','04','05')] == "./.") < 5,]

    #### iterate through
    # SV_table - table of SVs
    # mummer - mapping table
    SV_table$ref <- NA
    SV_table$ref_pos <- NA

    for (scf in unique(SV_table$scf)) {
        # scf <- '4_Tte_b3v06_scaf000001'
        scf_sv <- SV_table[SV_table$scf == scf,]
        scf_aln <- mummer[mummer$q_name == scf,]

        refs <- c()
        poss <- c()

        for ( line in 1:nrow(scf_sv) ) {
            dist_to_alignments <- abs(scf_aln$q_start - scf_sv$pos[line])
            if ( min(dist_to_alignments) < treshold_for_chaining){
                window_aln <- scf_aln[which.min(dist_to_alignments), ]
                refs <- c(refs, window_aln$r_name)
                poss <- c(poss, window_aln$r_start + (window_aln$q_start - scf_sv$pos[line]))
            } else {
                refs <- c(refs, NA)
                poss <- c(poss, NA)
            }
        }

        SV_table[SV_table$scf == scf, 'ref'] <- refs
        SV_table[SV_table$scf == scf, 'ref_pos'] <- poss
    }

    # mean(is.na(SV_table$ref))
    # I placed ~75% of variants

    SV_table$lg <- gsub("^([[:alnum:]]*).*", "\\1", SV_table$ref)

    # remove unassigned variants
    SV_table <- SV_table[!is.na(SV_table$lg),]
    table(SV_table$lg)
    # separate them by LG
    LGs <- paste0("lg", c(1:13, "X"))
    SV_table_per_LG <- lapply(LGs, function(lg){ SV_table[SV_table$lg == lg,] } )

    ####
    reference$lg <- gsub("^([[:alnum:]]*).*", "\\1", reference$scf)

    rownames(reference_map) <- reference_map$scf


    for ( i in 1:12 ){
        lg <- LGs[i]
        LG_ref <- reference[reference$lg == lg, ]

        lg_len <- sum(LG_ref$len)
        LG_homo_SVs <- rep(0, lg_len / 1e6)
        LG_hetero_SVs <- rep(0, lg_len / 1e6)
        LG_both_SVs <-  rep(0, lg_len / 1e6)

        SVs <- SV_table_per_LG[[i]]
        homoz <- rowSums(SVs[,c('00','01','02','03','04','05')] == "1/1")
        het <- rowSums(SVs[,c('00','01','02','03','04','05')] == "0/1")
        SVs$homo <- homoz > 0 & het == 0
        SVs$het <- het > 0 & homoz == 0

        for (line in 1:nrow(SVs)){
            adj_pos_from <- round((SVs$ref_pos[line] + reference_map[SVs$ref[line], 'adj']) / 1e6)
            adj_pos_to <- round((SVs$ref_pos[line] + (SVs$pos[line] - SVs$end[line]) + reference_map[SVs$ref[line], 'adj']) / 1e6)
            if ( SVs$homo[line] ){
                LG_homo_SVs[adj_pos_from:adj_pos_to ] <- LG_homo_SVs[adj_pos_from:adj_pos_to] + 1
            }
            if ( SVs$het[line] ){
                LG_hetero_SVs[adj_pos_from:adj_pos_to ] <- LG_hetero_SVs[adj_pos_from:adj_pos_to] + 1
            }
            if ( SVs$homo[line] == SVs$het[line] ){
                LG_both_SVs[adj_pos_from:adj_pos_to ] <- LG_hetero_SVs[adj_pos_from:adj_pos_to] + 1
            }
        }

        file <- paste0('figures/SVs/', sp, '_', lg, '.pdf')

        pdf(file, width = 8, height = 4)
            ylim = range(c(LG_homo_SVs, LG_hetero_SVs, LG_both_SVs), na.rm = T)
            plot(LG_homo_SVs, ylim = ylim, type = 'l', lty = 1, main = paste(timemas$names[sp_i], lg) , xlab = 'Position on chromosome [ Mbp ]', ylab = '# found SVs', col = 1)
            lines(LG_hetero_SVs, lty = 2, col = 2)
            lines(LG_both_SVs, lty = 3, col = 4)
            legend('topright', bty = 'n', c('homozygous in all', 'heterozygous in all', 'both'), lty = 1:3, col = c(1,2,4))
        dev.off()
    }
}


