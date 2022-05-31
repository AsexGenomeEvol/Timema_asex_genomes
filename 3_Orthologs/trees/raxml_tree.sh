## tree

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### strat.
## use 10sp orths (genome-only and all) to draw a tree
## Use keiths alignments (mcoffee) - these are back-translated form the aa so are codon aligned
## cat aligns --> filter? --> RAxML
## q - do with an outgroup??


############################################################################################
### Final Alignment is here 3_Orthologs/trees/tree_data/align_five_pair_HOGs_all.fa 
### get aligmnents of orthologs from selectome (1-to-1, five pair, codons)
## join alignments
## filter alignment with gblocks

~/Gblocks_0.91b/Gblocks align_five_pair_HOGs_all.fa -t=c -e=-gb1 -b4=12 
#Original alignment: 10503174 positions
#Gblocks alignment:  2377398 positions (22 %) in 7788 selected block(s)


############################################################################################
## draw tree
module add Bioinformatics/Software/vital-it
module load Phylogeny/raxml/8.2.12

cat > codon_specify_all.txt ## 3_Orthologs/trees/tree_data/codon_specify_all.txt
DNA, codon1 = 1-2377398\3
DNA, codon2 = 2-2377398\3
DNA, codon3 = 3-2377398\3

raxmlHPC -m GTRGAMMA -n five_pair_HOGs_all_1000boot -s align_five_pair_HOGs_all_gb.fa -p 12345  -# 1000 -q codon_specify_all.txt -f a -x 12345 -c 40 -T 1


### Tree here:
# 3_Orthologs/trees/output/RAxML_bipartitions.five_pair_HOGs_all_1000boot





