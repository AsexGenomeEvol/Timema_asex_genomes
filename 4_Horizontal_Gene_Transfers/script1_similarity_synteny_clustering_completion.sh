#!/bin/bash     

################## SCRIPT STANDARDIZED, FOR ARTHROPODS (cf 'outline_pipeline.pdf' for an overview of the pipeline steps & parameters)  ############################

# this script :
#			1 - performs for each Timema species (independently) the similarity and synteny steps -> HGT candidates
#			2 - clustering of all HGT candidates detected in the 10 Timema species -> putative HGT families
#			3 - completion of all HGT families ('rescue' of homologs not yet identified in the genomic scaffolds) -> completed putative HGT families
#			4 - retrieve the reference ('blast-hitting') sequences for each HGT family -> ready for alignment

# (the subsequent alignments and tree reconstruction have to be run afterwards (see "script2_phylogenetic_validation.sh"))


# script run on the 10 Timema species.
# this version: parallelized on 40 CPU (can be modified).
# this script should be launched from the folder in which all results will be written / output.

# run with the b3v04 genome assemblies ( = BEFORE contamination removal with blobtools, to be sure not to loose any HGT-containing scaffold)
# b3v04 assemblies = deposited NCBI assemblies + sequences removed from the NCBI assemblies (available here: https://github.com/AsexGenomeEvol/Timema_asex_genomes/tree/main/4_Horizontal_Gene_Transfers/contamination_sequences).



########### SOFTWARE REQUIREMENTS:
# required softwares: DIAMOND blast / Gmap / Silix / samtools / MAFFT / faidx (pyfaidx) / blast+
# required scripts (in ~/scripts/): NonOverlap_Merger_v2.pl / identity.r (courtesy of Paul Simion)


########### DATA REQUIREMENTS:

#		- a genome assembly (scaffolds) for each species
#		- the corresponding set of predicted coding sequences (CDS) (both as cDNA and protein sequences) (+ corresponding .fai files)
#		   WARNING: headers should be consistent between the CDS_cDNA and CDS_protein fasta files (at least the first field of the header should be exactly similar in both files).

#		- a protein reference database (for the DIAMOND blast step)(should be DIAMOND-formatted, i.e. "database.dmnd"). This database needs to be balanced and to represent the actual diversity of living organisms.
#		   (in this analysis, I used a custom database: 'customDB_2017_09_arthropod.dmnd'. For details, see Francois et al. 2020, G3, https://doi.org/10.1534/g3.119.400758 )

#		- a table 'species' listing the info for species to be analyzed (ex: Tbi sex T. bartmani) (no header) ( !! the 1st field should correspond to the abbreviation used in genomic files !!!)



######## OUTPUT:
# 'summ_after_gmap1' -> summary file of the number of candidate HGT for each species
# folder fasta_for_align/ -> HGT-families (completed, fasta files including 50 or 100 ref sequences, ready for alignment)
# log file of the script



######## FOLDERS TO BE DEFINED BY THE USER: [ !! you need to adapt this to your architecture / data !! ]

# specify the working rep where all output will be written:
rep="$(pwd)"
# specify where are the required files [table 'species' & reference database .dmnd] for all HGT analyses:
hgt="$(echo ~/timema/hgt)"
# path to the fasta files with CDS sequences (as cDNA and protein = files rna & pep respectively + their fai files):
rna="$(echo ~/timema/transcriptomes)"
# path to the genomic scaffolds (assembly WITH contaminants left)
dna="$(echo ~/timema/genomes/final_v3/with_contam)"


######## WARNING: in this script, files are called this way to match my data structure:
# genome="$(ls "$dna"/*_"$sp"_b3v04.fa)"
# cds="$(ls "$rna"/"$sp"_alliso_v2_minRKPM_2__300_longest_iso.fa.transdecoder.cds.longest)"
# pep="$(echo "$rna"/"$sp"_alliso_v2_minRKPM_2__300_longest_iso.fa.transdecoder.pep.longest)"
# reference database for  DIAMOND blast: "/home/clementine/timema/hgt/customDB_2017_09_arthropod" (called line 136)
# also: line 702 there is a command line to adapt : need the file name for PEP !
########## you need to adapt these to your data  !!!!!!!! (attention: they are used several times !!)





######## save the date of analyses and the software versions:
DATE=`date +%Y-%m-%d`
echo "Analyses were run the $DATE" > settings
echo `diamond --version` >> settings
echo `gmap --version | sed '3q;d'` >> settings
echo `silix --version` >> settings
echo `samtools --version | head -n 2` >> settings
pathblast="$(ls ~/software/ncbi-blast*/bin/blastp | cut -f1-6 -d '/')"
export PATH=$PATH:$pathblast
echo `blastp -version` >> settings
echo "MAFFT version is:" >> settings
mafft --version 2>>settings
echo "faidx (pyfaidx) version is:" >> settings
echo `faidx --version` >> settings





######## prepare the summary files which will be filled species by species during the analyses:

# 'summ_after_gmap1' shows the resulst after the blast step (similarity) ND after the Gmap step (synteny) = the most comprehensive table:
printf ""species"\t"sex"\t"N_orf"\t"N_orf_with_tax"\t"N_hgt1"\t"N_hgt1_eubact"\t"N_hgt1_arch"\t"N_hgt1_protist"\t"N_hgt1_fungi"\t"N_hgt1_plant"\t"N_hgt1_uncertainonlynm"\t"N_orf_arthropod"\t"N_orf_metazoa"\t"N_orf_mapped"\t"percent_orf_mapped"\t"N_hgt2"\t"N_hgt_eubact2"\t"N_hgt_arch2"\t"N_hgt_protist2"\t"N_hgt_fungi2"\t"N_hgt_plant2"\t"N_hgt_uncertainonlynm2"\n" > "$rep"/summ_after_gmap1

categ="archaea eubacteria fungi plant protist uncertainonlynm"
categ_all="archaea eubacteria fungi plant protist uncertainonlynm arthropod othermetazoa"
categ_HGT="archaea eubacteria fungi plant protist"


##################################### FIRST loop on all species:

for sp in `cut -f1 "$hgt"/species`
do

echo "------> Start processing $sp."

long="$(grep "$sp" "$hgt"/species | cut -f3)"
sex="$(grep "$sp" "$hgt"/species | cut -f2)"

if [ ! -d "$sp" ]; then
mkdir "$sp"
fi
cd "$sp"/

genome="$(ls "$dna"/*_"$sp"_b3v04.fa)"
cds="$(ls "$rna"/"$sp"_alliso_v2_minRKPM_2__300_longest_iso.fa.transdecoder.cds.longest)"
pep="$(echo "$rna"/"$sp"_alliso_v2_minRKPM_2__300_longest_iso.fa.transdecoder.pep.longest)"
norf="$(grep -c '>' "$pep")"


# warning if the size filter went wrong (in the transcriptome):
min="$(cut -f1-2 "$pep".fai |  sort -k2,2n | head -n 1 | cut -f2)"
if (( $(echo "$min < 99" | bc -l) )); then
echo "PROBLEM: minimum size of transcripts < 99 amino acids"
fi



######## 1. BLAST step:

echo "Start the first BLAST step for $sp."

# first create the DIAMOND alignment archive (DAA) output (stores all information):
# option --seg (yes,no) Enable SEG masking of low complexity segments in the query
diamond blastp -d /home/clementine/timema/hgt/customDB_2017_09_arthropod -q "$pep" --threads 40 --seg yes --evalue 1E-5 --max-target-seqs 150 --more-sensitive -o "$sp"_diamond_blastp.daa --outfmt 100 1>"$sp"_diamond_blastp.log 2>"$sp"_diamond_blastp.err

# then create, from the DAA, the accurate tabular output: (here, we require 14 fields):
diamond view --daa "$sp"_diamond_blastp.daa --out "$sp"_diamond_blastp.tab --threads 40 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen sseq 1>"$sp"_diamond_view.log 2>"$sp"_diamond_view.err
# this outformat corresponds to the 12 default + 2 additional fields:
# qlen 		Query sequence length
# sseq 		Aligned part of subject sequence (here the gaps in the alignment are not indicated)


blast="$(echo "$sp"_diamond_blastp.tab)"

mkdir "$sp"_blast

# first filter: minimum of 40% identity, min 75AA alignment length and Evalue <= E-10 -> list of all ORFs with at least 1 confident blast hit:
LC_ALL=C sort -k 11,11 -g $blast | LC_ALL=C awk '$11<=1e-10{print}' | LC_ALL=C sort -k 3,3 -g | LC_ALL=C awk '$3>=40{print}' | LC_ALL=C sort -k 4,4 -g | LC_ALL=C awk '$4>=75{print}'| cut -f1 | sort | uniq > "$sp"_blast/"$sp"_list_withhit
n_with_blast_hit="$(wc -l < "$sp"_blast/"$sp"_list_withhit)"

comm -3 <(grep '>' "$pep" | sed 's/>//g' | cut -d ' ' -f1 | sort) <(sort "$sp"_blast/"$sp"_list_withhit) > "$sp"_blast/"$sp"_nohit
n_no_blast_hit="$(wc -l < "$sp"_blast/"$sp"_nohit)"


# prepare the files for each taxonomic group (see 'categ_all') and remove them if they already exist:
#  2>/dev/null redirects the error message (if the file does not exist) to the 'dump' folder = avoid polluting the screen.
for cat in $categ_all
do
rm "$sp"_blast/"$sp"_"$cat" 2>/dev/null; touch "$sp"_blast/"$sp"_"$cat"
done
#categ_all="archaea eubacteria fungi plant protist uncertainonlynm arthropod othermetazoa"



####### parse the BLAST output, based on defined criteria for taxonomic assignment:
# on the 10 best unique hits (unique -> keep only 1 hit per reference sequence, as to not give too much weight to a single seq which might be erroneous).
# hits from at least 2 species are required (otherwise: discard this ORF)
# we choose to be stringent for the 'arthropod' genes, because later used to confirm that 'dubious' genes are actually integrated IN the host genome]
# arthropod = only metazoa hits and at least 50% arthropod (including at least 2 arthropod species)
# othermetazoa (actually: these are not really  used afterwards) = only metazoa hits and less than 50% arthropod (including at least 2 metazoa species)
# for other groups: criteria are relaxed to account for potential contamination (from symbionts, diet, or other) in genome assemblies.
# plant / eubacteria / archaea / fungi / protist = at least 70% of the corresponding taxonomic group
# uncertain = all ORFs meeting none of the above requirements for a clear taxonomic assignment
# (including) uncertainonlynm = same as above BUT only non-metazoa hits.


# loop on all ORFs with at least 1 confident blast hit ( with 1e-10 evalue, 40% identity and alignment length of 75 AA at least):


for contig in `cat "$sp"_blast/"$sp"_list_withhit`
do

# 1. create 'contigtmp' with the 10 best unique hits that we will be using:		
# awk '!x[$0]++' enables to remove the redundant ref sequences while keeping the order (by evalue) = keeping the best hit of this ref sequence
grep $contig $blast | LC_ALL=C sort -k 3,3 -g | LC_ALL=C awk '$3>=40{print}' | LC_ALL=C sort -k 4,4 -g | LC_ALL=C awk '$4>=75{print}'| LC_ALL=C sort -k 11,11 -g | LC_ALL=C awk '$11<=1e-10{print}' | cut -f2 | awk '!x[$0]++' | head -n 10 > contigtmp
# list of the 'detailed group' of the 10 hits (ie arthropod / eubacteria etc)
cut -d '|' -f2 contigtmp > temp
# uniq ( metazoa / nonmetazoa ) of the 10 hits:
cut -d '|' -f1 contigtmp | sort | uniq > temp2
# cut the 'detailed | species' field  (ex: arthropod|drosophila_melanogaster)
cut -d '|' -f2,3 contigtmp | sort | uniq > temp3

occ=`cat temp | wc -l`
fifty=`awk "BEGIN {print "$occ"/2}"`
seventy=`awk "BEGIN {print "$occ"*7/10}"`

# set all 'counters' at 0 [ in the 10 best hits: arth = nber of arthropod HITS / narth = number of arthropod SPECIES ]
arth="$(echo "0")"; narth="$(echo "0")"

# number of different species in the 10 best hits:
nspec="$(cut -d '|' -f2 temp3 | sort | uniq | wc -l)"

# CASE: arthropod + othermetazoa
if  grep -Fxq "arthropod" temp; then
arth=`sort temp | uniq -c | grep 'arthropod' | awk '{print $1}'`
fi

narth="$(grep -w 'arthropod' temp3 | cut -d '|' -f2 | sort | uniq | wc -l)"

if  [[ `uniq temp` = "arthropod" ]] && (( $(echo "$narth >= 2" | bc -l) )); then
echo "$contig" >>"$sp"_blast/"$sp"_arthropod
fi
if [[ `cat temp2` = "metazoa" ]] && (( $(echo "$arth >= $fifty" | bc -l) )) && grep -Fxq "othermetazoa" temp; then
echo "$contig" >> "$sp"_blast/"$sp"_arthropod
fi
if [[ `cat temp2` = "metazoa" ]] && (( $(echo "$arth < $fifty" | bc -l) )) && grep -Fxq "othermetazoa" temp && (( $(echo "$nspec >= 2" | bc -l) )); then
echo "$contig" >> "$sp"_blast/"$sp"_othermetazoa
fi


## CASE = suspected HGT : loop on the different 'categ_HGT' of HGT ( = archaea / eubacteria / fungi / plant / protist ) :

for cat in $categ_HGT
#categ_HGT="archaea eubacteria fungi plant protist"
do
# set the counters at 0 for each category:
doubt="$(echo "0")"; ndoubt="$(echo "0")"

if  grep -Fxq "$cat" temp; then

doubt=`sort temp | uniq -c | grep -w "$cat" | awk '{print $1}'`
ndoubt="$(grep -w "$cat" temp3 | cut -d '|' -f2 | sort | uniq | wc -l)"
if  [[ `uniq temp` = "$cat" ]] && (( $(echo "$ndoubt >= 2" | bc -l) )); then
echo "$contig" >>"$sp"_blast/"$sp"_"$cat"
fi
if (( $(echo "$doubt >= $seventy" | bc -l) )) && (( $(echo "$ndoubt >= 2" | bc -l) )) && [[ `uniq temp` != "$cat" ]]; then
echo "$contig" >> "$sp"_blast/"$sp"_"$cat"
fi

fi

done


done



# create the list of ORFs with no 'reliable' taxonomic assignment = "uncertain"
cat "$sp"_blast/"$sp"_arthropod "$sp"_blast/"$sp"_othermetazoa > "$sp"_blast/"$sp"_all_tax_assigned
for cat in $categ_HGT
do
	cat "$sp"_blast/"$sp"_"$cat" >> "$sp"_blast/"$sp"_all_tax_assigned
done
comm -3 <(sort "$sp"_blast/"$sp"_list_withhit) <(sort "$sp"_blast/"$sp"_all_tax_assigned) > "$sp"_blast/"$sp"_uncertain


# then, create the sub-list of these 'unreliable-taxonomy' ORFs which have ONLY non-metazoa hits = "uncertainonlynm" => these will be kept !
for contig in `cat "$sp"_blast/"$sp"_uncertain`
do
nspec="$(echo "0")"
grep $contig $blast | LC_ALL=C sort -k 3,3 -g | LC_ALL=C awk '$3>=40{print}' | LC_ALL=C sort -k 4,4 -g | LC_ALL=C awk '$4>=75{print}'| LC_ALL=C sort -k 11,11 -g | LC_ALL=C awk '$11<=1e-10{print}'| cut -f2 | awk '!x[$0]++' | head -n 10 > contigtmp
cut -d '|' -f1 contigtmp | sort | uniq > temp2
cut -d '|' -f2,3 contigtmp | sort | uniq > temp3
nspec="$(cut -d '|' -f2 temp3 | sort | uniq | wc -l)"
if [[ `cat temp2` = "nonmetazoa" ]] && (( $(echo "$nspec >= 2" | bc -l) )); then
echo "$contig" >> "$sp"_blast/"$sp"_uncertainonlynm
fi
done



# cleaning the temp files:
rm contigtmp 2>/dev/null; rm temp 2>/dev/null; rm temp2 2>/dev/null; rm temp3 2>/dev/null


# also generate a summary file "$sp"_all_tax with the taxonomic assignment for all ORFs (when it was achievable):
rm "$sp"_all_tax 2>/dev/null; touch "$sp"_all_tax
# for the taxonomic groups: archaea / eubacteria / fungi / plant / protist / uncertainonlynm / arthropod / othermetazoa [uncertain is removed]
#categ_all="archaea eubacteria fungi plant protist uncertainonlynm arthropod othermetazoa"
for cat in $categ_all
do
for orf in `cat "$sp"_blast/"$sp"_"$cat"`
do
printf "$orf\t$cat\n" >> "$sp"_all_tax
done
done
# then, remove all 'uncertain', because tricky to manipulate afterwards ! (doublons with 'uncertainonlynm')
sed -i '/uncertain$/d' "$sp"_all_tax



# summarize these intermediary results:
# (these will be printed on the joint summary file 'summ_after_gmap1' with blast & Gmap intermediary results.)
norfwithtax="$(wc -l < "$sp"_blast/"$sp"_all_tax_assigned)"
norfarthropod="$(wc -l < "$sp"_blast/"$sp"_arthropod)"
o="$(wc -l < "$sp"_blast/"$sp"_othermetazoa)"
norfmetazoa=`awk "BEGIN {print "$norfarthropod"+"$o"}"`
nhgteub="$(wc -l < "$sp"_blast/"$sp"_eubacteria)"
nhgtarch="$(wc -l < "$sp"_blast/"$sp"_archaea)"
nhgtprot="$(wc -l < "$sp"_blast/"$sp"_protist)"
nhgtfung="$(wc -l < "$sp"_blast/"$sp"_fungi)"
nhgtplant="$(wc -l < "$sp"_blast/"$sp"_plant)"
nhgt=`awk "BEGIN {print "$nhgteub"+"$nhgtarch"+"$nhgtprot"+"$nhgtfung"+"$nhgtplant"}"`
nhgtunconlynm="$(wc -l < "$sp"_blast/"$sp"_uncertainonlynm)"


echo "End of the BLAST step for $sp."


################ 2. Gmap step (synteny):

echo "Start the Gmap1 step for $sp."


# Rq: here, we are still in the "$sp" folder
# reminder: rna="$(echo ~/timema/transcriptomes)" and dna="$(echo ~/timema/genomes/final_v3/with_contam)"
# Reminder: how we set the file name for the genomic scaffolds:
#genome="$(ls "$dna"/*_"$sp"_b3v04.fa)"
#cds="$(ls "$rna"/"$sp"_alliso_v2_minRKPM_2__300_longest_iso.fa.transdecoder.cds.longest)"

mkdir "$sp"_gmap1



#### 2a. for Gmap, first need to 'index' the reference genome (if it does not exist yet), for each species (here, one genomic scaffold per FASTA entry):
# redirect the standard output / error to "$sp"_refgenome/gmap_build.log and gmap_build.err
# this 'indexed' ref genome is stored in "$dna"/
if [ ! -d "$dna"/"$sp"_refgenome ]; then
echo "Start indexing the reference genome $sp for Gmap"
cd "$dna"; mkdir "$sp"_refgenome; cd "$sp"_refgenome
gmap_build -D "$dna"/"$sp"_refgenome -d "$sp" "$genome" 1>"$dna"/"$sp"_refgenome/gmap_build.log 2>"$dna"/"$sp"_refgenome/gmap_build.err
cd "$rep"/"$sp"
fi



#### 2b. map all ORFs (for which taxonomic assignment was achievable / reliable in the previous blast step) to the ref genomic scaffolds:
# we have defined a minimum length of 100 bp and min identity of 95% for this alignment step.


# create the list of to-be-Gmapped ORFs (for which taxonomic assignment was achievable in the previous blast step), if it doesn't exist yet:
if [ ! -f "$sp"_gmap1/"$sp"_orf_with_tax ]; then
cat "$sp"_blast/"$sp"_all_tax_assigned "$sp"_blast/"$sp"_uncertainonlynm > "$sp"_gmap1/"$sp"_orf_with_tax
fi

# retrieve all the corresponding cDNA sequences (=> 'all_tax_assigned' + 'uncertainonlynm') in a fasta file:
cat "$sp"_gmap1/"$sp"_orf_with_tax | xargs -n 1 samtools faidx "$cds" > "$sp"_cds


# Gmap: map all cDNA on the 'indexed' genomic scaffolds (min 95% identity).
# option -n 0 to allow for (what they call) 'chimeras' (= 1 ORF spanning 2 or more scaffolds) while not allowing for multimapping (a portion of 1 ORF map to several genomic locations)
# -n = Maximum number of paths to show (default 5).If you want a single alignment plus chimeric alignments, then set this to be 0.
# the rationale here is that a given transcript is supposed to originate from only one location in the genome.
# and the problem is that we observed that multimapping creates some inconsistent links between scaffolds of incongruent taxonomy ....
# We redirect the SAM std output to "$sp"_gmap.sam and std error to "$sp"_gmap.err
gmap -D "$dna"/"$sp"_refgenome -d "$sp" -B 5 -t 40 -n 0 -f samse --no-sam-headers --min-identity=0.95 "$sp"_cds 1>"$sp"_gmap1/"$sp"_gmap.sam 2>"$sp"_gmap1/"$sp"_gmap.err

# cleaning:
rm "$sp"_cds


# Reminder: to check / visualize the alignment (in SAM format):use sam2pairwise to see in a human-readable view the alignment (open the output with geany)
#grep 'Tsi_TRINITY_DN29652_c0_g3_i2|m.6970' ex_sam | ~/software/sam2pairwise-master/src/sam2pairwise > ~/tmp; geany ~/tmp





#### 2c. create a table summarizing the info for each mapped cDNA:  ORF id / scaffold on whit it mapped / taxonomy of the ORF
# when no mapping was achievable for a given ORF: indicated in the SAM file by $2==4
# require a min. length of 100 bp (of the aligned segments (~ exons)) (~ exon length in metazoa & empirical on some visual tests)

rm "$sp"_gmap1/"$sp"_nomapping 2>/dev/null; touch "$sp"_gmap1/"$sp"_nomapping
rm "$sp"_gmap1/"$sp"_mapping 2>/dev/null; touch "$sp"_gmap1/"$sp"_mapping


for orf in `cat "$sp"_gmap1/"$sp"_orf_with_tax`
do
tax="$(echo NA)"
grep -w "$orf" "$sp"_gmap1/"$sp"_gmap.sam > samtemp
# if the ORF was not reliably mapped:
if  [[ `cut -f2 samtemp` == 4 ]]; then
echo "$orf" >> "$sp"_gmap1/"$sp"_nomapping
else
# assign the taxonomic group to each orf based on the blast step ("$sp"_all_tax)
# here the taxonomic groups are: archaea / eubacteria / fungi / plant / protist / arthropod / othermetazoa / uncertainonlynm
tax="$(grep -w "$orf" "$sp"_all_tax | cut -f2)"


## now we summarize this mapping/tax info (for each ORF), in the file "$sp"_mapping:
# here, applies to the cases where the ORF mapped on 1 or more scaffolds:
nmap="$(cat samtemp | wc -l)"
for i in $(seq 1 1 $nmap)
do
# retrieve the focal line:
sed "${i}q;d" samtemp > subtemp
# then, retrieve the MD string to calculate the number of aligned bases (matches=md1 + mismatches=md2):
# attention! the SAM output of Gmap is weird: when a given cDNA maps to > 1 scaffold, the MD tag is in the 13th field instead of the 12th...
md="$(tr '[ \t]' '[\n\n]' < subtemp | grep '^MD' | cut -d':' -f3)"
# md1 = matches = sum of the numbers indicative of matches in the MD string
md1="$(echo $md | grep -o -E '[0-9]+'| awk '{ SUM += $1} END { print SUM }')"
# md2 = mismatches = sum of the letters indicative of mismatches in the MD string
md2="$(echo $md | sed 's/[0-9]//g' | wc -c)"
# require a min length of alignment (for the mapping parts ~ exons) of 100 bp:
if (( $(echo "$md1 + $md2 - 1 >= 100" | bc -l) )); then
scaff="$(cut -f3 subtemp)"
# print one line per scaffold on which this ORF mapped: ORF id / scaffold ID / taxonomy assigned to this ORF
printf "$orf\t$scaff\t$tax\n" >> "$sp"_gmap1/"$sp"_mapping
fi
rm subtemp
done
fi
rm samtemp
done


# loop to check if there is any 'NA' taxonomy remaining => problem
if grep -Fwxq "NA" <(cut -f3 "$sp"_gmap1/"$sp"_mapping); then
echo "There is problem of taxonomic assignment for some ORF in the species $sp"
fi





#### 2d. sort the HGT-candidates based on this table:
# keep only those which are located on 'arthropod-like' scaffolds (= on which at least 1 'arthropod-vertically-descended' ORF mapped)


# create the list of HGT-candidates (after the blast step, including uncertainonlynm), if it doesn't exist yet:
rm "$sp"_gmap1/"$sp"_candid_after_blast 2>/dev/null; touch "$sp"_gmap1/"$sp"_candid_after_blast
#categ="archaea eubacteria fungi plant protist uncertainonlynm"
for cat in $categ
do
cat "$sp"_blast/"$sp"_"$cat" >> "$sp"_gmap1/"$sp"_candid_after_blast
done

# prepare the files for each category (archaea / eubacteria / fungi / plant / protist / uncertainonlynm) (remove them if they already exist):
for cat in $categ
do
rm "$sp"_gmap1/"$sp"_"$cat" 2>/dev/null; touch "$sp"_gmap1/"$sp"_"$cat"
done



### scann all HGT-candidates, decide if they are to be kept or not, and assign them to their taxonomic group (eubacteria / plant / ...):
for cand in `cat "$sp"_gmap1/"$sp"_candid_after_blast`
do
# retrieve the taxonomic group of this HGT-candidate:
taxcand="$(grep -w "$cand" "$sp"_all_tax | cut -f2)"
# on which scaffold(s) did it map (can be more than 1):
scaff2="$(grep -w $cand "$sp"_gmap1/"$sp"_mapping | cut -f2)"
# then retrieve the taxonomic assignment of ALL ORFs which have been mapped to this (or these) given scaffold(s):
cat "$sp"_gmap1/"$sp"_mapping | grep -w "$scaff2" | cut -f3 | sort | uniq > temp
## now, we filter those HGT-candidates which shared at least 1 scaffold with another 'arthropod-like' gene: [ contamination filter, sensitive to N50 !! ]
if (grep -Fwxq "arthropod" temp); then
echo "$cand" >> "$sp"_gmap1/"$sp"_"$taxcand"
fi
rm temp
done


# create the list of all HGT-candidates (including 'uncertainonlynm') surviving this Gmap filter: "$sp"_candid_after_gmap:
rm "$sp"_gmap1/"$sp"_candid_after_gmap 2>/dev/null; touch "$sp"_gmap1/"$sp"_candid_after_gmap
for cat in $categ
do
cat "$sp"_gmap1/"$sp"_"$cat" >> "$sp"_gmap1/"$sp"_candid_after_gmap
done


# summarize these intermediary results and print them in the summary file, jointly with those of the Blast step:
n_orf_mapped="$(cut -f1 "$sp"_gmap1/"$sp"_mapping | sort | uniq | wc -l)"
n_orf_to_be_mapped="$(wc -l < "$sp"_gmap1/"$sp"_orf_with_tax)"
p_orf_mapped="$(echo $(( n_orf_mapped *100 / n_orf_to_be_mapped )))"
nhgteub2="$(wc -l < "$sp"_gmap1/"$sp"_eubacteria)"
nhgtarch2="$(wc -l < "$sp"_gmap1/"$sp"_archaea)"
nhgtprot2="$(wc -l < "$sp"_gmap1/"$sp"_protist)"
nhgtfung2="$(wc -l < "$sp"_gmap1/"$sp"_fungi)"
nhgtplant2="$(wc -l < "$sp"_gmap1/"$sp"_plant)"
nhgt2=`awk "BEGIN {print "$nhgteub2"+"$nhgtarch2"+"$nhgtprot2"+"$nhgtfung2"+"$nhgtplant2"}"`
nhgtunconlynm2="$(wc -l < "$sp"_gmap1/"$sp"_uncertainonlynm)"


printf "$sp\t$sex\t$norf\t$norfwithtax\t$nhgt\t$nhgteub\t$nhgtarch\t$nhgtprot\t$nhgtfung\t$nhgtplant\t$nhgtunconlynm\t$norfarthropod\t$norfmetazoa\t$n_orf_mapped\t$p_orf_mapped\t$nhgt2\t$nhgteub2\t$nhgtarch2\t$nhgtprot2\t$nhgtfung2\t$nhgtplant2\t$nhgtunconlynm2\n" >> "$rep"/summ_after_gmap1

#printf ""species"\t"sex"\t"N_orf"\t"N_orf_with_tax"\t"N_hgt1"\t"N_hgt1_eubact"\t"N_hgt1_arch"\t"N_hgt1_protist"\t"N_hgt1_fungi"\t"N_hgt1_plant"\t"N_hgt1_uncertainonlynm"\t"N_orf_arthropod"\t"N_orf_metazoa"\t"N_orf_mapped"\t"percent_orf_mapped"\t"N_hgt2"\t"N_hgt_eubact2"\t"N_hgt_arch2"\t"N_hgt_protist2"\t"N_hgt_fungi2"\t"N_hgt_plant2"\t"N_hgt_uncertainonlynm2"\n" > "$rep"/summ_after_gmap1



echo "End of the Gmap step for $sp, switch to the next species"

cd "$rep"


done   ### end of the loop on the species !!!!



##################################### Generate some summary tables:

echo "End of the first loop on all species, generating some summary files of all species"

# concatenate all blast output tables (for all species) in a summary table 'concat_all_blast' covering all species:
cd "$rep"
rm concat_all_blast 2>/dev/null; touch concat_all_blast

for sp in `cut -f1 "$hgt"/species`
do
	cat "$sp"/"$sp"_diamond_blastp.tab  >> concat_all_blast
done

# there is a bug with some lines containing only 'blosum62'... (not sure why, likely introduced by diamond view). Just remove these lines, if they exist:
sed -i '/^blosum62/d' concat_all_blast




##################################### Clustering of all HGT-candidates in families with Silix:

echo "Start the clustering of HGT-candidates in families with Silix."
cd "$rep"


# first, need a fasta file with the AA sequences of all HGT-candidates from all species: 'seqs_aa_all_candid_1.fa':
rm seqs_aa_all_candid_1.fa 2>/dev/null; touch seqs_aa_all_candid_1.fa
for sp in `cut -f1 "$hgt"/species`
do
pep="$(ls "$rna"/"$sp"_alliso_v2_minRKPM_2__300_longest_iso.fa.transdecoder.pep.longest)"
cat "$sp"/"$sp"_gmap1/"$sp"_candid_after_gmap | xargs -n 1 samtools faidx "$pep" >> seqs_aa_all_candid_1.fa
done

n_candid1="$(grep -c '>' seqs_aa_all_candid_1.fa)"


mkdir silix; cd silix
# then, perform an all-against-all BLAST (makeblastdb, then blastp, output in tabular format):
mkdir blast
makeblastdb -in "$rep"/seqs_aa_all_candid_1.fa -dbtype prot -out blast/blastdb -title blastdb 1>blast/makedb.log
blastp -num_threads 40 -query "$rep"/seqs_aa_all_candid_1.fa -db blast/blastdb -evalue=1e-04 -outfmt=6 -out blast/allvsallblastp 2>blast/blastp.err


# finally, clustering with Silix (all default parameters, except for min identity of 85% which is quite stringent ! ):
# the rationale here is to 'overcluster' the families, because we prefer to oversplit them than to merge some families which shouldn't have been merged.
# indeed, what we look for are 'mono-species' families, so if we don't find any while oversplitting, it's highly unlikely that there is any !
silix -n -i 0.85 -r 0.80 -l 100 -m 0.50 "$rep"/seqs_aa_all_candid_1.fa blast/allvsallblastp -f FAM > silix_results 2>silix.err

# retrieving family sizes:
pathsilix="$(ls ~/software/silix*/utils/silix-fsize | cut -f1-6 -d '/')"
"$pathsilix"/silix-fsize silix_results > silix_fam_size
cut -f2 silix_fam_size | sort -n | uniq -c > distrib_fam1


rm temp_sp 2>/dev/null; touch temp_sp
for fam in `cut -f1 silix_results | sort | uniq`
do
grep -w $fam silix_results | cut -f2 | cut -f1 -d '_' | sort | uniq | wc -l >> temp_sp
done
sort temp_sp | uniq -c > distrib_fam1_sp
rm temp_sp


n_fam="$(wc -l < silix_fam_size)"

# Splitting sequences in multiple fasta files corresponding each to one FAM (some clusters have a size of 1):
mkdir fam_silix
silix-split -o fam_silix -n 1 -p "seqs" "$rep"/seqs_aa_all_candid_1.fa silix_results 2>silix-split.err


cd "$rep"

echo "Silix clustered $n_candid1 HGT-candidates into $n_fam families. See their distribution (n families / n HGT-sequences within each fam):"
cat silix/distrib_fam1
echo "and their distribution in terms of species (n families / n species within each fam):"
cat silix/distrib_fam1_sp


echo "End of the first clustering step."





##################################### Loop on the families for their 'completion' / rescue of some candidates:

echo "Start the completion of the Silix families."

# for each HGT-fam, try to complete it by searching for this HGT-candidate in the genomes of ALL species (with Gmap)
# Reminder: all FAM are (in fasta, pep seqs) in the folder "$rep"/silix/fam_silix/
# here, we create a new folder to store the FASTA files for all completed families: "$rep"/silix/fam_silix_completed/
# all completion analyses / details are stored in the folders "$rep"/silix/completion/famXX/"

## Rq: this step controls for variations in N50 among species:
#-> one HGT shared between spA (bad N50) & spB (great N50) will likely be kept in B / lost in A by the synteny filter (gmap1), but will be rescued in A here.
# ( ! except for 'private HGT' in low-N50 species which are hopelessly lost...)


# **********************************************************************
##### OVERVIEW OF THE COMPLETION PROCESS:

##### for a given family F, try to complete it with some seqs from the species S, when this species S being searched is NOT YET in this family F:
 
# if more than 1 seq in the fam F, select the longest seq => our 'reference-HGT-candidate' for this fam F (when the spS was already in the famF, we select the longest seq of the spS in famF).
# Gmap this 'reference-HGT-candidate' for this fam F (cDNA of the ORF) on the genomic scaffolds of the species S. (min %ident 85 % for the aligned parts ( ~ exons))
# retrieve the proteic sequences of all matching exons, concatenated on each scaffold (options -Q -Y) -> 1 concatenated AA sequence per 'matching' scaffold (! same ID before '|' for all !)
# exclude all "concatenated AA sequences" shorter than 33 AA (~ the threshold of 100 nt used previsouly) -> min of 33 AA aligned / scaffold.
# with MAFFT (--auto), align these 'rescued' AA sequences onto the corresponding AA seq of the 'reference-HGT-candidate' for this fam F.
# needed by the following script: seq in one line / remove the ref seq from the alignment / create a file 'notification' / all the considered concatenated AA sequence have the same ID before the first '|'
# use the script 'NonOverlap_Merger_v2.pl' (P. Simion) to detect if 'significant' overlap (> 7 AA) between these rescued sequences & retrieve the focal seqs based on that (different decision rules):

# IF NO OVERLAP (or no more than 25 AA) between 'rescued':
	# we assume this is a unique gene (most plausible hyp): the script concatenates all "rescued" AA seqs in the order given by the alignment => create a 'consensus' sequence (containing '-' that are kept)
	# notification "fusion" -> add this consensus concatenated seq as 'spS_rescued_famF' to fam F.

# IF OVERLAP between 'rescued' (at least on 26 AA): 
	# we assume there are several paralogs: notification "overlap" -> need to perform an additional step.
	# (in bash) scan the alignment to find the max number of rescued sequences overlapping at a given position = N = the number of paralogs in species S for this fam.
	# (if there are several overlaps of N sequences, select the longest).
	# for each such paralog spanning this 'best' overlap -> retrieve the concatenated-exon sequence and add as 'spS_rescued_N_famF' to fam F (N = number)
	# (do not try to merge with other rescued fragments outside of this overlapping region, bcz impossible to know which fragment goes with another => would create chimeras)


#### When the species S being searched is ALREADY in the family F: we need to search for possible other duplicates, BUT take care not to rescue the gene already present in the fam F = would artefactually 'duplicate' it !
# this turns out to be quite difficult! We consider 2 cases:
# - when there is only 1 seq of spS in famF -> we can device something (because we can confidently identify the 'reference seq' (here, the longest seq of the spS in famF) in the rescued seqs -> 100 %identity):
#		if only 1 sequence rescued (no overlap): do nothing as we are rescuing the same (single) ref seq that was used to search into the genome.
#		if there are N >1 rescued seqs (overlap): out of these N seqs, remove the seq with the highest %identity (relative to the ref seq => most likely the same) and rescue the (N-1) others.
# - when there are > 1 seq of spS in famF -> give up, not possible / highly tricky to confidently identify the "primary" seqs (already present in the famF) in the rescued seqs.


######### reminder: here what matters the most is to rescue at least 1 seq per species, to know if a given HGT event/family is species-specific, or shared. Given the parameters we chose, we MAY miss some duplicates, but this is not our focus.
# **********************************************************************


cd silix/
mkdir fam_silix_completed
mkdir completion

cut -f1 silix_fam_size > list_fam		# create the list of the Silix families to be completed

## start of the LOOP ON THE FAMILIES
for fam in `cat list_fam`
do

cp "$rep"/silix/fam_silix/seqs_"$fam".fa fam_silix_completed/seqs_"$fam".fasta

cd completion
mkdir "$fam"_completion; cd "$fam"_completion


# first, retrieve the cDNA sequences of the HGT-candidates (ORF) of this family:
grep -w "$fam" "$rep"/silix/silix_results | cut -f2 > listorf_"$fam"
rm seqs_"$fam" 2>/dev/null; touch seqs_"$fam"

for sp in `cut -f1 "$hgt"/species`  
do
# reminder:
#genome="$(ls "$dna"/*_"$sp"_b3v04.fa)"
transcriptome="$(ls "$rna"/"$sp"*cds.longest)"
cat listorf_"$fam" | grep "$sp" | xargs -n 1 samtools faidx $transcriptome >> seqs_"$fam"
done
# retrieve the species name of the candidates primarily present in this family:
cut -f1 -d '_' listorf_"$fam"  > listsp_"$fam"
nsp="$(grep -c "$sp" listsp_"$fam")"

## start of the LOOP ON THE SPECIES
# in each family F, try to complete it with some seqs from each of the 10  species:

for sp in `cut -f1 "$hgt"/species`
do
samtools faidx seqs_"$fam"
## select which sequence in the famF will be used as a reference:
# if the species S being searched was NOT YET in this famF: select the longest seq as our 'reference-HGT-candidate' for this fam F
if grep -vq $sp listsp_"$fam"; then
idref="$(sort -k2 -n seqs_"$fam".fai | tail -n1 | cut -f1)"
samtools faidx seqs_"$fam" $idref > seqs_"$fam"_ref
else	# otherwise:  the species S being searched was ALREADY present in this famF: select the longest seq of this spS in famF as our 'reference-HGT-candidate' for this fam F
idref="$(grep "$sp" seqs_"$fam".fai | sort -k2 -n | tail -n1 | cut -f1)"
samtools faidx seqs_"$fam" $idref > seqs_"$fam"_ref
fi
# Gmap this 'reference-HGT-candidate' for this fam F (cDNA of the ORF) on the genomic scaffolds of the species S:
# min %ident 85 % for the aligned parts ( ~ exons)
# here, multimapping is allowed as we look also for duplicates within the genomes.
# retrieve the proteic sequences of all matching exons, concatenated on each scaffold  -> 1 concatenated AA sequence per 'matching' scaffold (! same ID before '|' for all !)
# (here, options  -P to print protein sequence and  -Y to translate cDNA with corrections for frameshifts (equivalent to MACSE?)).
# reminder: the 'indexed' ref genomes are stored in "$dna"/sp_refgenome
gmap -D "$dna"/"$sp"_refgenome -d "$sp" -B 5 -t 40 -Q -Y --min-identity=0.85 seqs_"$fam"_ref > gmap_"$sp".out 2>"$sp"_gmap.err
# remove all "concatenated AA sequences" shorter than 33 AA (~ corresponding to the threshold of 100 nt used previosuly):
if  [[ -s gmap_"$sp".out ]]; then
awk 'BEGIN{RS=">";FS="\n"}NR>1{seq="";for (i=2;i<=NF;i++) seq=seq""$i; print ">"$1"\n"seq}' gmap_"$sp".out | awk '!/^>/ {next}{ getline seq} length(seq) >= 33 {print $0 "\n" seq}' > tmp.fa
rm gmap_"$sp".out; mv tmp.fa gmap_"$sp".out
fi

# keep on only if Gmap found something (the output file 'gmap_"$sp".out' is not empty):
if  [[ -s gmap_"$sp".out ]]; then
# need to modify the headers of the output protein seqs, because Gmap gives the (same) ID to all, those of the 'reference-HGT-candidate': add "_rescuedN" to the end of the header:
awk '/^>/{gsub(/$/,"_rescued"i++" ");}1' gmap_"$sp".out > gmap_"$sp".out2
# add the PROTEIN sequence of the 'reference-HGT-candidate':
spref="$(echo "$idref" | cut -f1 -d '_')"
samtools faidx $rna/"$spref"_alliso_v2_minRKPM_2__300_longest_iso.fa.transdecoder.pep.longest "$idref" >> gmap_"$sp".out2
# with MAFFT (--auto), align these 'rescued' AA sequences onto the corresponding AA seq of the 'reference-HGT-candidate' for this fam F.
mafft --auto gmap_"$sp".out2 > mafft_"$sp".aln 2>"$sp"_mafft.err
# required for the subsequent script: sequences on one single line / remove the ref seq from the alignment fasta file / create 2 files 'notification' & 'identity' (info about the overlap and %id resp.)
faidx mafft_"$sp".aln -g "rescued" | awk 'BEGIN{RS=">";FS="\n"}NR>1{seq="";for (i=2;i<=NF;i++) seq=seq""$i; print ">"$1"\n"seq}' > mafft_"$sp".aln2
rm notification_"$sp" 2>/dev/null; touch notification_"$sp"
rm identity_"$sp" 2>/dev/null; touch identity_"$sp"
# use the script 'NonOverlap_Merger_v2.pl' (P. Simion) to detect if 'significant' overlap (> 25 AA) between these rescued sequences & retrieve the focal seqs based on that (different decision rules):
perl ~/scripts/NonOverlap_Merger_v2.pl mafft_"$sp".aln2 Merger_"$sp".out notification_"$sp" identity_"$sp" 1>Merger_"$sp".log

## now, different actions depending if "fusion" (= no overlap) or "overlap" (at least on 26 AA):

## IF NO OVERLAP (or no more than 25 AA) between 'rescued': [notification_sp = "fusion"] -> we assume this is a unique gene (most plausible hyp)
# the script has concatenated all "rescued" AA seqs in the order given by the alignment => create a 'consensus' sequence (containing '-' that are kept)
if grep -Fxq "fusion" notification_"$sp"; then
#  when the species S being searched was NOT YET in this famF: the consensus concatenated seq is added as 'spS_rescued' to the fasta file of the famF in the folder fam_silix_completed/
if ! grep -q $sp listsp_"$fam"; then
sed "s/>.*$/>$sp\_rescued/g" Merger_"$sp".out | cat >> "$rep"/silix/fam_silix_completed/seqs_"$fam".fasta
fi

#  when the species S being searched was ALREADY present in this family F: do nothing, as we just retrieved the same (single) sequence already present in the famF.
fi			# end of the loop for the notification "fusion"


## IF OVERLAP  between 'rescued' (at least on 26 AA):  [notification_sp = "overlap"] -> we assume there are several paralogs.
if grep -Fxq "overlap" notification_"$sp"; then
# scan the alignment to find the max number of rescued sequences overlapping at a given position = N = the ESTIMATED number of paralogs in spS for this famF (can be underestimated).
grep '>' mafft_"$sp".aln2 | sed 's/>//g' > list_rescued_"$sp"
samtools faidx mafft_"$sp".aln2
lalign="$(cut -f2 mafft_"$sp".aln2.fai | sort | uniq | tail -n 1)"		# length of this alignment
rm overlap 2>/dev/null; touch overlap	# file where we store the number of overlapping rescued seqs (one line per position in the alignment)
for i in $(seq 1 1 $lalign)			# loop on all positions of the alignment
do
rm temp 2>/dev/null; touch temp
for id in `cat list_rescued_"$sp"`
do
samtools faidx mafft_"$sp".aln2 $id:$i-$i | tail -n1 >> temp
done
sed 's/-//g' temp | wc -w >> overlap		# retrieve the number of sequences overlapping at this position ('-' indicate a missing position for this rescued seq / a letter indicate an AA)
rm temp
done
maxoverlap="$(sort -n overlap | uniq | tail -n 1)"	# N = maximum number of overlapping sequences for this alignment
# if there are several overlaps of N rescued seqs, select the longest; and retrieve the (final) position (in the alignment) of this longest overlap:
grep -n $maxoverlap overlap | sed 's/:[0-9]*//g' > overlap2	# get all positions in the alignment where the number of overlapping sequences = N
noverlap2="$(wc -l < overlap2)"
tmp="$(echo 0)"; winner="$(echo 0)"
for  i in $(seq 2 1 $noverlap2)
do
j=$(($i - 1))
ival="$(sed ''$i'q;d' overlap2)"
jval="$(sed ''$j'q;d' overlap2)"
jvalplus=$((jval + 1))
if [ "$ival" == "$jvalplus" ]; then	# if the value of the line i equal the value of the line j + 1 (= sequential suite of nbers) -> we are WITHIN an overlap: increase the length of this overlap by 1 and keep on
tmp=$((tmp+1))
else						# if not: we are moving from an overlap to another: 
if [[ "$tmp" -gt "$winner" ]]; then	# if this overlap (that we are leaving) is longer than the last one (or 0 if it's the first): save this max overlap length (as $winner) and the (finale) position of this max overlap
winner="$(echo $tmp)"
position="$(echo "$jval")"
fi
tmp="$(echo 0)"					# reset the counter (length of the overlap) to 0, before starting the new overlap
fi
done
# final check: if there was only 1 overlap, or if the alignment ends on the longest overlap, the above loop would fail to get the accurate numbers: so check:
if [[ "$tmp" -gt "$winner" ]]; then
position="$(echo "$jval")"
fi

#  when the species S being searched was NOT YET in this famF:
# for all N paralog spanning this 'best' overlap -> retrieve the concatenated-exon sequence and add them as 'spS_rescued_N' to the fasta file of the famF in the folder fam_silix_completed/ (N = number)
if ! grep -q "$sp" listsp_"$fam"; then
rm rescued 2>/dev/null; touch rescued	# file where we store the ID + fasta seq of the 'chosen' overlapping seqs that we will rescue
for id in `cat list_rescued_"$sp"`
do
letter="$(samtools faidx mafft_"$sp".aln2 $id:$position-$position | tail -n1)"
if grep -q [A-Z] <<< "$letter" ; then
samtools faidx mafft_"$sp".aln2 $id >> rescued
fi
done
sed "s/>.*$/>$sp/g" rescued | awk '/^>/{gsub(/$/,"_rescued"i++" ");}1' | cat >> "$rep"/silix/fam_silix_completed/seqs_"$fam".fasta
rm rescued
fi
#  when the species S being searched was ALREADY present in this famF (only 1 seq of this spS)
# if there were more than 1 seq of spS in famF: give up.
# if there was only 1 seq of spS in famF: out of these N seqs, remove the seq with the highest %identity (relative to the ref seq => most likely the same) and rescue the (N-1) others.
# => add them as 'spS_rescued_i' to the fasta file of the famF in the folder fam_silix_completed/
if grep -q "$sp" listsp_"$fam" && [[ $nsp == 1 ]];then
rm rescued 2>/dev/null; touch rescued		# file where we store the ID of the N overlapping seqs
for id in `cat list_rescued_"$sp"`
do
letter="$(samtools faidx mafft_"$sp".aln2 $id:$position-$position | tail -n1)"
if grep -q [A-Z] <<< "$letter" ; then
echo "$id" >> rescued		# get the ID of the N overlapping seqs we are rescuing. (one has to be removed)
fi
done
cp mafft_"$sp".aln aligntmp

# calculate (in R) the pairwise %id of all sequences (including the ref): call the script 'identity.r' to do the job:
~/scripts/identity.r

tail -n+2 distances > distances_"$sp"
# remove the non-overlapping seqs -> in the remaining N rescued overlapping sequences: identify the seq with the highest %id, remove it and get the (N-1) ID of the seqs to be actually rescued:
grep -f rescued distances_"$sp" | sort -n -k2 -t ' ' | head -n -1 | cut -f1 -d ' ' > idrescued
# then, retrieve the (N-1) corresponding protein seqs and add them as 'spS_rescued_i' to the fasta file of the famF in the folder fam_silix_completed/
grep -f rescued distances_"$sp" | sort -n -k2 -t ' ' | head -n -1 | cut -f1 -d ' ' | cat |  xargs -n 1 samtools faidx mafft_"$sp".aln2 > rescued2
sed "s/>.*$/>$sp/g" rescued2 | awk '/^>/{gsub(/$/,"_rescued"i++" ");}1' | cat >> "$rep"/silix/fam_silix_completed/seqs_"$fam".fasta
rm rescued; rm rescued2; rm aligntmp; rm distances
fi
rm overlap; rm overlap2

fi		# end of the loop for the notification "overlap"

fi

done     # end of the LOOP ON THE SPECIES

cd "$rep"/silix

done # end of the LOOP ON THE FAMILIES


## calculate some stats:
rm distrib_fam_completed 2>/dev/null; touch distrib_fam_completed
for file in `ls "$rep"/silix/fam_silix_completed/seqs*fasta`
do
n="$(grep -c '>' $file)"
echo $n >> distrib_fam_completed
done
n_candid2="$(awk '{ sum += $1 } END { print sum }' distrib_fam_completed)"
n_rescued=`awk "BEGIN {print "$n_candid2"-"$n_candid1"}"`

sort -n distrib_fam_completed | uniq -c > distrib_fam2


rm distrib_fam_completed_sp 2>/dev/null; touch distrib_fam_completed_sp
for file in `ls "$rep"/silix/fam_silix_completed/seqs*fasta`
do
n="$(grep '>' $file | cut -f1 -d '_' | sort | uniq | wc -l)"
echo $n >> distrib_fam_completed_sp
done

sort -n distrib_fam_completed_sp | uniq -c > distrib_fam2_sp


rm fam_completed_sp_detail 2>/dev/null; touch fam_completed_sp_detail
for file in `ls "$rep"/silix/fam_silix_completed/seqs*fasta`
do
grep '>' $file | cut -f1 -d '_' | sed 's/>//g' >> distrib_fam_completed
done


echo "AFTER completion of the $n_fam families: $n_candid2 HGT-sequences (this step rescued $n_rescued candidate-sequences)."
echo "See their distribution (n families / n HGT-sequences within each fam):"
cat distrib_fam2
echo "and their distribution in terms of species (n families / n species within each fam):"
cat distrib_fam2_sp
echo "Detail of the number of HGT-derived sequences for each species after completion:"
sort distrib_fam_completed | uniq -c

cd "$rep"




##################################### Loop on the 'final' families: get the  50 best 'blast-hitting' reference seqs (for subsequent automated alignment)

echo "Start to retrieve the 'blast-hitting' reference seqs for each family."
echo "(= build the FASTA files containing the AA sequences of Timema + ref species from customDB -> alignments for phylogenetic validation)."

cd "$rep"
# we retrieve max. 50 ref sequences (from customDB) per family:
mkdir fasta_for_align; cd fasta_for_align
mkdir ref50

### rationale (when there are n > 1 timema seqs in the family):
# we retrieve the 50 best hits for each of the n timema seqs independently
# merge the n * 50 best | sort | uniq (use awk '!_[$2]++' to uniq on the 2nd column while keeping the good order which conveys the info on the Evalue of each hit)
# and then take the 50 best ( -> to avoid redundancy !!!)

#### !!!! retrieve ONLY the 'matching/aligned' part of the customDB sequences blasting on each HGT-candidate !!!
# previous comparisons (whole vs. truncated) on some families  seemed to indicate that alignments were better when truncated.

# minimum of alignment length at 50 AA (in diamond-blast, this 'length' corresponds to the number of aligned positions in the query = excludes all gaps)

for fam in `cat "$rep"/silix/list_fam`		# loop on the families
do
file="$(echo seqs_"$fam")"
cp "$rep"/silix/fam_silix_completed/"$file".fasta ref50/"$fam".fasta
grep '>' "$rep"/silix/fam_silix/"$file".fa | sed 's/>//g' > tmp		# IDs of the timema seqs of this family (not the 'rescued' seqs = no blast result for them)

rm contigtmp 2>/dev/null; touch contigtmp
for orf in `cat tmp`
do
# get the 50 UNIQUE 'best' hits for each timema seq of this famF (min. 40% identity / E-10 evalue / 75 AA length)
# (uniq while keeping the order ~ Evalue: awk '!_[$2]++' )
# concatenate these best hits for all timema seqs in famF:
grep -w "$orf" "$rep"/concat_all_blast | LC_ALL=C sort -k 3,3 -g | LC_ALL=C awk '$3>=40{print}' | LC_ALL=C sort -k 4,4 -g | LC_ALL=C awk '$4>=75{print}' | LC_ALL=C sort -k 11,11 -g | LC_ALL=C awk '$11<=1e-10{print}'  | awk '!_[$2]++' | head -n 100 >> contigtmp
done

LC_ALL=C sort -k 11,11 -g contigtmp | awk '!_[$2]++' | head -n 50 > contigtmp50

while read LINE; do
echo $LINE | cut -f2 -d ' ' >> ref50/"$fam".fasta
echo $LINE | cut -f14 -d ' ' >> ref50/"$fam".fasta
done < contigtmp50


sed -i 's/^nonmetazoa/>nonmetazoa/g' ref50/"$fam".fasta
sed -i 's/^metazoa/>metazoa/g' ref50/"$fam".fasta

rm contigtmp; rm contigtmp50; rm tmp
done

cd "$rep"

echo "End of this script: families of HGT-candidates have been identified, completed and the corresponding FASTA files are available at "$rep"/fasta_for_align/"


#### END of the pipeline => HGT-candidates families (fasta files) ready for alignment then phylogenetic validation !!!



