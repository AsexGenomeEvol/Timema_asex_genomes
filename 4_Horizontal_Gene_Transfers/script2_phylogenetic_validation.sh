#!/bin/bash     

################## script for phylogenetic validation of candidate HGT families (cf 'outline_pipeline.pdf' for an overview of the pipeline steps & parameters)  ############################

# this script starts from the completed putative HGT families (identified with the script "script1_similarity_synteny_clustering_completion.sh").
#			1 - automated alignment (with the 50 best blast hits with the reference database)
#			2 - Tree reconstruction 

## Softwares: MAFFT (--version v7.313), HMMcleaner, RAXML
## scripts GapCleaner2, HMMcleanAA.pl & last_step3.py: courtesy of Emeric Figuet



############### Alignments   #####################################################################################################################

## Loop on the HGT families :
for famlong in /home/clementine/Timema/new_analyses_completion_aout2019/fasta_for_align/FAM*.fasta
do

fam=`echo $famlong | rev | cut -f1 -d '/' | rev`
prefix=`echo $fam | cut -f1 -d '.' | rev | cut -f1 -d '/' | rev`

# align with MAFFT:
mafft --auto "$famlong" > "$fam"_aln
# shorten the headers to keep only ID (HMMcleaner only considers the Nth first characters, so it could wrongly consider 2 seqs to have identical ID)
sed -i 's/>[a-z]*|[a-z]*|[a-zA-Z0-9_-;.: ]*|/>/g' "$fam"_aln

# HMMCleaner -> mask 'weak' positions in the alignment:
# --del-char : character ('@') to replace bad quality sites (will be replaced by '-' afterwards).
# 12 : stringency parameter (default = 10, I increased it a bit to be less stringent, because we are dealing with highly divergent seqs).
~/scripts/HMMcleanAA.pl --del-char @ "$fam"_aln 12
sed -i 's/@/-/g' "$prefix"_Hmm12.fasta

# remove the 'columns'/sites for which >50% of missing info
python3.5 ~/scripts/GapCleaner2 "$prefix"_Hmm12.fasta 0.5

# remove the sequences too short (< 25% of the alignment length after previous step)
~/scripts/last_step3.py $prefix

# finally: restore the full ID of the sequences(ID with 4 fields):
for seq in `grep '>' "$prefix".final.fasta | sed 's/>//g'`
do
temp=`grep "$seq" "$famlong"`
sed -i "s/>$seq/$temp/g" $prefix.final.fasta
done

done





###############  Phylogenetic reconstruction   #################################################################################################


mkdir -p intermediaires_trees final_trees

## Loop on the HGT families :
for fam in FAM*.final.fasta
do
prefix=`basename "$fam" .final.fasta`

sed "s/:/_/g" "$fam" | sed "s/;/_/g" > "$prefix"_temp.fa

# builds the 'best' tree: -> RAxML_bestTree.prefix
raxmlHPC-PTHREADS-SSE3 -m PROTGAMMALGX -s "$prefix"_temp.fa -n $prefix -p 5 -# 50 -T 40
# separate run for the bootstraps : (option -b: to add bootstrap) -> RAxML_bootstrap.prefix_bs
raxmlHPC-PTHREADS-SSE3 -m PROTGAMMALGX -s "$prefix"_temp.fa -n ${prefix}_bs -p 5 -# 100 -b 123 -T 40
# To map the results of the bootstrap on the tree (besttree) -> RAxML_bipartitions.prefix.tr
raxmlHPC-PTHREADS-SSE3 -f b -z RAxML_bootstrap.${prefix}_bs -t RAxML_bestTree.$prefix -m PROTGAMMALGX -s "$prefix"_temp.fa -n $prefix.tr -T 40 

# Formatting (final files):
mv RAxML_bipartitions.$prefix.tr RAxML.${prefix}_FINAL.tr; mv RAxML.${prefix}_FINAL.tr final_trees/
mv ${prefix}*reduced RAxML_bootstrap.${prefix}_bs RAxML_*.$prefix.RUN.* RAxML_bipartitionsBranchLabels.$prefix.tr RAxML_info.$prefix* RAxML_bestTree.$prefix intermediaires_trees

rm "$prefix"_temp.fa

done

# Options:
#-#50 for the exploration (alternative starting tree): 30 or 50 is enough (number of alternative runs on distinct starting trees)
#-p 5: random seed
#-n: output file
#-b 123: seed for bootstrap
#-T: CPU







