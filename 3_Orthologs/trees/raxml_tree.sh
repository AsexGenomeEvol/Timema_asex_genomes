## tree

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### strat.
## use 10sp orths (genome-only and all) to draw a tree
## Use keiths alignments (mcoffee) - these are back-translated form the aa so are codon aligned
## cat aligns --> filter? --> RAxML
## q - do with an outgroup??


############################################################################################
### Final Alignment is here 3_Orthologs/trees/tree_data/align_five_pair_HOGs_genomeonly_gb.fa 
### get aligmnents of orthologs from selectome (1-to-1, five pair, genome copies only, codons)
## join alignments
## filter alignment with gblocks

~/Gblocks_0.91b/Gblocks align_five_pair_HOGs_genomeonly.fa -t=c -e=-gb1 -b4=12 
#Original alignment: 6272457 positions
#Gblocks alignment:  1523166 positions (24 %) in 4520 selected block(s)

############################################################################################
## draw tree
module add Bioinformatics/Software/vital-it
module load Phylogeny/raxml/8.2.12

cat > codon_specify_genome_only.txt ## 3_Orthologs/trees/tree_data/codon_specify_genome_only.txt
DNA, codon1 = 1-1523166\3
DNA, codon2 = 2-1523166\3
DNA, codon3 = 3-1523166\3

raxmlHPC -m GTRGAMMA -n five_pair_HOGs_genomeonly_1000boot -s align_five_pair_HOGs_genomeonly_gb.fa -p 12345  -# 1000 -q codon_specify_genome_only.txt -f a -x 12345 -c 40 -T 1


### Tree here:
# 3_Orthologs/trees/output/RAxML_bipartitions.five_pair_HOGs_genomeonly_1000boot





