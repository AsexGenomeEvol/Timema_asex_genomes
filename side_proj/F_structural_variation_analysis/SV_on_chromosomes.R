library(AsexStats)
library(RColorBrewer)

# load species genome info
scf_len_files <- paste0("data/", timemas$codes, "/reference/", timemas$codes, "_b3v08_scf.lengths")
coord1_files <- paste0('data/', timemas$codes, '/genome_aln/to_ref/',
                       timemas$codes, '_b3v08_to_ref.1coords')


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
# sp_i <- 7
source('../timema_assembly/scripts/prepare_reference.R')

get_SV_on_LGs_table <- function(sp_i){
    # parameters
    treshold_for_chaining <- 10000
    cov_treshold <- 0.6

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

    SV_table$lg <- gsub("^([[:alnum:]]*).*", "\\1", SV_table$ref)
    print(mean(is.na(SV_table$lg)))
    # I placed ~75% of variants of Tce

    # remove unassigned variants
    SV_table <- SV_table[!is.na(SV_table$lg),]
    # table(SV_table$lg)

    homoz <- rowSums(SV_table[,c('00','01','02','03','04','05')] == "1/1")
    het <- rowSums(SV_table[,c('00','01','02','03','04','05')] == "0/1")
    SV_table$homo <- homoz > 0 & het == 0
    SV_table$het <- het > 0 & homoz == 0
    SV_table
}

SV_tables <- lapply(1:10, get_SV_on_LGs_table)
# SAVE IT!!!
for (i in 1:10) {
    tab_to_save <- paste0("stats/", timemas$code[i], "_SVs_on_LGs.tsv")
    write.table(SV_tables[[i]], tab_to_save, quote = F, row.names = F, sep = '\t')
}

# SV_table_per_LG <- lapply(LGs, function(lg){ SV_table[SV_table$lg == lg,] } )

LGs <- paste0("lg", c(1:13, "X")) # a vector of LGs
reference$lg <- gsub("^([[:alnum:]]*).*", "\\1", reference$scf)
lg_sizes <- sapply(unique(reference$lg), function(x){ sum(reference[reference$lg == x, 'len']) } )
lg_adj <- c(0,cumsum(lg_sizes)[1:12])
rownames(reference_map) <- reference_map$scf
reference_map$lg <- gsub("^([[:alnum:]]*).*", "\\1", reference_map$scf)
reference_map$mahatten_adj <- reference_map$adj + lg_adj[factor(reference_map$lg, LGs)]
genome_len <- lg_adj[13]

resolution <- 5e6

# # By Myles Harrison
# addalpha <- function(colors, alpha=1.0) {
#   r <- col2rgb(colors, alpha=T)
#   # Apply alpha
#   r[4,] <- alpha*255
#   r <- r/255.0
#   return(rgb(r[1,], r[2,], r[3,], r[4,]))
# }

plot_SVs_on_LS <- function(sp_i, exclude_zero = T, resolution = 1e6, lines = F){
    SV_table <- SV_tables[[sp_i]]
    sp <- timemas$names[sp_i]
    SV_table <- SV_table[SV_table$lg != "lgX",]

    genome_homo_SVs <- rep(0, genome_len / resolution)
    genome_hetero_SVs <- rep(0, genome_len / resolution)
    # SV_table
    for (line in 1:nrow(SV_table)){
        adj_pos_from <- round((SV_table$ref_pos[line] + reference_map[SV_table$ref[line], 'mahatten_adj']) / resolution)
        adj_pos_to <- round((SV_table$ref_pos[line] + (SV_table$pos[line] - SV_table$end[line]) + reference_map[SV_table$ref[line], 'mahatten_adj']) / resolution)
        if ( SV_table$homo[line] ){
            genome_homo_SVs[adj_pos_from:adj_pos_to ] <- genome_homo_SVs[adj_pos_from:adj_pos_to] + 1
        }
        if ( SV_table$het[line] ){
            genome_hetero_SVs[adj_pos_from:adj_pos_to ] <- genome_hetero_SVs[adj_pos_from:adj_pos_to] + 1
        }
    }

    # BrBG
    if ( sp_i %% 2 == 1){
        pal <- brewer.pal(5, "YlGnBu")[c(3,5)]
    } else {
        pal <- brewer.pal(5, "YlOrRd")[c(3,5)]
    }
    # pal <- addalpha(pal, alpha)


    ylim = range(c(genome_homo_SVs, genome_hetero_SVs), na.rm = T)
    plot(NULL, xlim = c(1, genome_len / resolution), ylim = ylim, pch = 20,
         main = timemas$names[sp_i], xaxt = "n", bty = 'n', xlab = '', ylab = '', cex.axis = 1.4, cex.main = 1.6)
         # xlab = 'linage group [ Mbp ]', ylab = '# found SVs'
    xtick <- lg_adj / resolution
    axis(side = 1, at = xtick, labels = FALSE)
    text(x = (lg_adj[1:12] + lg_adj[2:13]) / (2 * resolution),  par("usr")[3],
         labels = LGs[1:12], pos = 1, xpd = TRUE, cex = 1.3)
    for (i in seq(1, 12, by = 2)) {
        rect(lg_adj[i] / resolution, -2, lg_adj[i + 1] / resolution, ylim[2], col = 'lightgrey', border = F)
    }

    if ( exclude_zero ){
        genome_homo_SVs[genome_homo_SVs == 0] <- NA
        genome_hetero_SVs[genome_hetero_SVs == 0] <- NA
    }

    if ( lines ){
        lines(genome_hetero_SVs, col = pal[2], lwd = 1.6)
        lines(genome_homo_SVs, col = pal[1], lwd = 1.6)
    } else {
        points(genome_hetero_SVs, pch = 20, col = pal[2])
        points(genome_homo_SVs, pch = 20, col = pal[1])
    }
        # legend('topright', bty = 'n', c('found only homozygous', 'found only heterozygous'), pch = 20, col = c(1,2))
}


for ( i in 1:10 ){
    sp <- timemas$codes[i]
    filename <- paste0('figures/SVs_mapped_to_LGs/', sp, '_manhatten_wo_zero.pdf')

    pdf(filename, width = 14, height = 4)
        plot_SVs_on_LS(i, T, 5e6, T)
    dev.off()
}

# mapped
# 1 - c(0.4411688, 0.4106708, 0.5809524, 0.4451531, 0.4884993, 0.6849438, 0.2665134, 0.3607531, 0.4903439, 0.962963)


# Tms colour morphs
# subset lg8
# SV_table <- SV_table[SV_table$lg == "lg8", ]
# subset only the region 38M - 52M
# SV_table_colour_locus <- SV_table[SV_table$ref_pos > 3.8e6 & SV_table$ref_pos < 5.2e6,]
# sapply(SV_table_colour_locus[,4:9], table)