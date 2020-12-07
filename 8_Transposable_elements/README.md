# Repeat identification and TE detection

This is how repeats were detected and annotated

## mining repeats from raw reads

Using dnaPipeTE v2.5 on serv04. dnaPipeTE had to be tweaked by Patrick to be able to run.


Raw read data put in scratch/local

- genome size must be specified
- parameters ok for most organizms: -genome_coverage 0.5 -sample_number 4

run on dee_serv04, because might consume up to 400Gb memory and highjack all CPUs


```
./DNApipeTE_run.sh

```

Output contains:

- annotated.fasta (contains all repeats that could be annotated as TEs)
- unannoted.fasta (contains all repeats that are not annotated)

For heterozygosity (theta) estimates, the two files were merged and the header reformatted to RepeatMasker standard and run with TE divergence (-div) 10, 20 ,30:

```
for f in *.fasta; do cat $f | awk '/^>/ $1 ~ "LTR" { printf("%s#%s/%s\n",$0,"LTR","LTR");next;} { print $0;}' | awk '/^>/ $1 ~ "LINE" { printf("%s#%s/%s\n",$0,"LINE","LINE");next;} { print $0;}' | awk '/^>/ $1 ~ "SINE" { printf("%s#%s/%s\n",$0,"SINE","SINE");next;} { print $0;}' | awk '/^>/ $1 ~ "DNA" { printf("%s#%s/%s\n",$0,"DNA","DNA");next;} { print $0;}' | awk '/^>/ $1 ~ "Helitron" { printf("%s#%s/%s\n",$0,"Helitron","Helitron");next;} { print $0;}' | awk '/^>/ $1 ~ "rRNA" { printf("%s#%s/%s\n",$0,"rRNA","rRNA");next;} { print $0;}' | awk '/^>/ $1 ~ "Low_complexity" { printf("%s#%s/%s\n",$0,"Low_complexity","Low_complexity");next;} { print $0;}' | awk '/^>/ $1 ~ "Satellite" { printf("%s#%s/%s\n",$0,"Satellite","Satellite");next;} { print $0;}' | awk '/^>/ $1 ~ "Simple_repeat" { printf("%s#%s/%s\n",$0,"Simple_repeat","Simple_repeat");next;} { print $0;}' | awk '/^>/ $1 ~ "na_comp" { printf("%s#%s/%s\n",$0,"Unknown","Unknown");next;} { print $0;}' > $f.rfmt; done

```

Two folders were made:
- raw_libraries (containing the repeat sequences)
- genomes (containing b3_v06 (or b3_v04) Timema genomes)

This script was run for RepeatMasker locally on dee_serv04 (check the parameters inside):

```
repeat_masking_local.sh
```

The output contains the soft-masked genome and a repeat overview table.


## mining repeats from assembles

additionally, RepeatModeler was run for repeat detection from assemblies from /scratch/local dee_serv04

```
repeatmodeler_run.sh
```

The output contains the file consensi.fa.classified
This file is a repeat library with annotated TEs.


## Clustering of repeat libraries

As both methods (raw and assembly based) of repeat detection have some shortcomings, we merged and clustered the two repeat libraries by 95% identity (the TE family threshold) using usearch

```
clustering.sh
```


## Annotation

The build-in annotation pipeline of dnaPipeTE and RepeatModeler are quite crude, which is why we will run PASTEC to classify TEs.

As there are many small sequences left after clustering, filtering > 500 (300) bp may be appropriate:
```
/scratch/beegfs/monthly/jbast/software/scripts/convert.py -i merged.TElib.fa.centroid95 -o merged.TElib.fa.centroid95min500 -s 500 -t fasta
```
Then, in preparation for pastec, header and sequences have to be reformatted to be interleaved and not contain special characters:
```
for f in *.centroid95min500; do python /scratch/beegfs/monthly/jbast/software/asv3.py -s rename_TE -i1 $f -o $f.cvt; done
```


If a database for host genes and/or rDNA is included, clean annotated host genes cleaned unknown // and get a TE sequences rDNA database (eg. SILVA):

```
#remove any TE-related annotations from these:
for f in *.fasta; do cat $f | grep -E "retro|transpos|transcriptase|unknown function|reverse|mobile element" | cut -d " " -f 1 | sed 's/>//g' > $f''.TEids; done
for f in *.fasta; do /scratch/beegfs/monthly/jbast/software/scripts/extract_contigs.py -i $f -o $f''.noTEs.fa -l $f''.TEids -r; done

#add a rDNA database (here SILVA used):
wget http://ftp.arb-silva.de/current/Exports/SILVA_132_SSURef_Nr99_tax_silva.fasta.gz     
wget http://ftp.arb-silva.de/current/Exports/SILVA_132_LSURef_tax_silva.fasta.gz          
cat SILVA_132_LSURef_tax_silva.fasta SILVA_132_SSURef_Nr99_tax_silva.fasta > rDNA_silva.fa  
```

PASTEClassifier: Please read readme_pastec. Specific folder structure has to be generated (with the CFG file in the same folder as samplenames and pastec_run.sh:

```
./pastec_run.sh
```
then, within each folder to run the complete classification pipe:

```
nohup PASTEClassifier_parallelized.py -i inpuFile.fa -C PASTEClassifier_parallelized.cfg -w -r  >& PASTEC.log &

```
This will submit many jobs to DEE-hugemem, so run it from /scratch/beegfs/monthls, not local !!

if something goes wrong, delete the MYSQL entries before restarting jobs:

```                                        
mysql --host=devbioinfo.unil.ch --user=synergia --password=synergia1234                
use synergia;                              
select * from jobs;                        

delete from jobs ; 
```                                        


Because PASTEClassifier did not annotate to TE family level, we blasted the repeat libraries classified by PASTEC
against the well curated Timema cristinae TE library generated for Soria-Carrasco et al. 2014.

```                                           
./blast_TEs_Tce.sh                               
```   

These blast hits have to be filtered according to TE classification standards:
filter ident%>80, alignment length >80, and the best best hit per contig was kept (merged) (bitscore and eval)
```
for f in *.blastn; do cat $f | awk '{if(($3>=80)&&($6>=80)){print($0)}}'| awk '{if(($3>=80)&&($6>=80)){print($0)}}' | sort -k1,1 -k10,10nr -k9,9n | sort -u -k1,1 --merge > $f''.filtered; done          
``` 

These blast annotations were added to the header of the PASTEc classification fasta files.
Needed some reformatting

add > in front of filtered blast                      
```
for f in *.filtered; do awk '{print(">"$0)}' $f > $f''2; done                                                 
```

add annotation to header if fits                      
```
for f in *.fa; do awk 'NR==FNR{a[$1]=$2;next}{print $0,"\t"a[$1]}' $f''.blastn.filtered2 $f > $f''.annot; done 
```

reformat annot headers to crop unnecessary parts      
```
for f in *.annot; do cat $f | sed 's/_.*\s/\t/' | sed 's/\s.*#/\t/' | sed 's/\s$//' > $f''2; done              
```

rename noCAT to Unknown if only classification (from PASTEC)
```
for f in *.annot2; do cat $f | sed 's/>noCat$/>Unknown/' > $f''.rn; done                                       
```

manually check all headers for both classifications and if they converge                  
```
cat *.annot2.rn | grep '>' | sort | awk '{if($0!=last)print($0);last=$0;}' > all.tes                           
```

convert to sed commands for all                             
```
cat all.tes | sed 's/\//\\\//' |sed 's/\s/\\s/' | sed "s/>/sed 's\//" | sed "s/$/$\/'/" > all.tes.seds 
```

Tweak this file to reannotate manually:
all.tes.seds

contains the manually selected and curated names with sed replacement command:
all.tes.cur

example:
sed 's/DTX-incomp\sDNA\/Kolobok$/DTX#DNA\/Kolobok/'


If conflict, PASTEC annotation was prioritized (with higher level of classification mostly)

Header naming was done with RepeatMasker standard naming but keeping the Wicker naming for elements:

Wicker#Repeatmasker e.g.:
DTA#DNA/hAT

```
all.tes.cur.sh
```

#the cleaned script is used like this:
```                                                           
./all.tes.cur.sh 1_Tdi_merged.TElib.fa.centroid95min500_negStrandReversed_WickerH.fa.annot2.rn 1_Tdi_merged.TElib.fa.centroid95min500.annot                                                                                                             
./all.tes.cur.sh 1_Tps_merged.TElib.fa.centroid95min500_negStrandReversed_WickerH.fa.annot2.rn 1_Tps_merged.TElib.fa.centroid95min500.annot                                                                                                             
./all.tes.cur.sh 2_Tcm_merged.TElib.fa.centroid95min500_negStrandReversed_WickerH.fa.annot2.rn 2_Tcm_merged.TElib.fa.centroid95min500.annot                                                                                                             
./all.tes.cur.sh 2_Tsi_merged.TElib.fa.centroid95min500_negStrandReversed_WickerH.fa.annot2.rn 2_Tsi_merged.TElib.fa.centroid95min500.annot                                                                                                             
./all.tes.cur.sh 3_Tce_merged.TElib.fa.centroid95min500_negStrandReversed_WickerH.fa.annot2.rn 3_Tce_merged.TElib.fa.centroid95min500.annot                                                                                                             
./all.tes.cur.sh 3_Tms_merged.TElib.fa.centroid95min500_negStrandReversed_WickerH.fa.annot2.rn 3_Tms_merged.TElib.fa.centroid95min500.annot                                                                                                             
./all.tes.cur.sh 4_Tbi_merged.TElib.fa.centroid95min500_negStrandReversed_WickerH.fa.annot2.rn 4_Tbi_merged.TElib.fa.centroid95min500.annot                                                                                                             
./all.tes.cur.sh 4_Tte_merged.TElib.fa.centroid95min500_negStrandReversed_WickerH.fa.annot2.rn 4_Tte_merged.TElib.fa.centroid95min500.annot                                                                                                            
./all.tes.cur.sh 5_Tge_merged.TElib.fa.centroid95min500_negStrandReversed_WickerH.fa.annot2.rn 5_Tge_merged.TElib.fa.centroid95min500.annot                                                                                                       
./all.tes.cur.sh 5_Tpa_merged.TElib.fa.centroid95min500_negStrandReversed_WickerH.fa.annot2.rn 5_Tpa_merged.TElib.fa.centroid95min500.annot
```


#sort fasta by contig header:                                 
sort_fasta.py                                                 
```
for f in *500.annot; do /scratch/beegfs/monthly/jbast/software/scripts/sort_fasta.py -i $f -o $f''.srt; done                
```

#then number and add species name                             
#as extra columns                                          
```
awk '{if(substr($1,1,1)==">"){if(lastname!=$1){i=1;}printf("%s\t%d\tTdi\n",$1,i++);lastname=$1;}else{print($0);}}'          
```

```
cat 1_Tdi_merged.TElib.fa.centroid95min500.annot.srt | awk '{if(substr($1,1,1)==">"){if(lastname!=$1){i=1;}printf("%s\t%d\tTdi\n",$1,i++);lastname=$1;}else{print($0);}}' > 1_Tdi_merged.TElib.fa.centroid95min500.annot.srt.nbr                        
cat 1_Tps_merged.TElib.fa.centroid95min500.annot.srt | awk '{if(substr($1,1,1)==">"){if(lastname!=$1){i=1;}printf("%s\t%d\tTps\n",$1,i++);lastname=$1;}else{print($0);}}' > 1_Tps_merged.TElib.fa.centroid95min500.annot.srt.nbr                        
cat 2_Tcm_merged.TElib.fa.centroid95min500.annot.srt | awk '{if(substr($1,1,1)==">"){if(lastname!=$1){i=1;}printf("%s\t%d\tTcm\n",$1,i++);lastname=$1;}else{print($0);}}' > 2_Tcm_merged.TElib.fa.centroid95min500.annot.srt.nbr                        
cat 2_Tsi_merged.TElib.fa.centroid95min500.annot.srt | awk '{if(substr($1,1,1)==">"){if(lastname!=$1){i=1;}printf("%s\t%d\tTsi\n",$1,i++);lastname=$1;}else{print($0);}}' > 2_Tsi_merged.TElib.fa.centroid95min500.annot.srt.nbr                        
cat 3_Tce_merged.TElib.fa.centroid95min500.annot.srt | awk '{if(substr($1,1,1)==">"){if(lastname!=$1){i=1;}printf("%s\t%d\tTce\n",$1,i++);lastname=$1;}else{print($0);}}' > 3_Tce_merged.TElib.fa.centroid95min500.annot.srt.nbr                        
cat 3_Tms_merged.TElib.fa.centroid95min500.annot.srt | awk '{if(substr($1,1,1)==">"){if(lastname!=$1){i=1;}printf("%s\t%d\tTms\n",$1,i++);lastname=$1;}else{print($0);}}' > 3_Tms_merged.TElib.fa.centroid95min500.annot.srt.nbr                        
cat 4_Tbi_merged.TElib.fa.centroid95min500.annot.srt | awk '{if(substr($1,1,1)==">"){if(lastname!=$1){i=1;}printf("%s\t%d\tTbi\n",$1,i++);lastname=$1;}else{print($0);}}' > 4_Tbi_merged.TElib.fa.centroid95min500.annot.srt.nbr                        
cat 4_Tte_merged.TElib.fa.centroid95min500.annot.srt | awk '{if(substr($1,1,1)==">"){if(lastname!=$1){i=1;}printf("%s\t%d\tTte\n",$1,i++);lastname=$1;}else{print($0);}}' > 4_Tte_merged.TElib.fa.centroid95min500.annot.srt.nbr                        
cat 5_Tge_merged.TElib.fa.centroid95min500.annot.srt | awk '{if(substr($1,1,1)==">"){if(lastname!=$1){i=1;}printf("%s\t%d\tTge\n",$1,i++);lastname=$1;}else{print($0);}}' > 5_Tge_merged.TElib.fa.centroid95min500.annot.srt.nbr                        
cat 5_Tpa_merged.TElib.fa.centroid95min500.annot.srt | awk '{if(substr($1,1,1)==">"){if(lastname!=$1){i=1;}printf("%s\t%d\tTpa\n",$1,i++);lastname=$1;}else{print($0);}}' > 5_Tpa_merged.TElib.fa.centroid95min500.annot.srt.nbr                        
```

#run RepeatMasker
```
repeatmasking_local.sh
```

#run repeat landscape scripts
```
repeat_landscapes.sh
```



#check if all TEs in TE librarary are found in the assemblies
```
blastn.sh
```
-> all there


#estimate TE content based on reads
Because not all repeats might be assembled, we base TE content estimates on the fraction of reads that map to TEs out of total mappable reads.
For this we use HTSeq-counts with the RepeatMasker and species specific TE library  generated gff annotation files and sorted bam files of bwa mapped reads to genomes.

#mapping
```
map_reads.sh
```

#gff files can be filtered, such that only TEs counted with speciefic minimum TE divergence and length
```
filter.sh
```

#HTSeq-counts
```
htseq.sh
```

#counts script to clean and generate R readable table
```
count.sh
```

#R script for count analyses
```
Timema_TEs.R
```

## RepeatMasker output to R plot (R dir)

1) Convert RepeatMasker output to table for R

```
python te.py -s repeatmasker2r -i1 <repeat_masker_file> -o <output>

python te.py -s repeatmasker2r -i1 1_Tps_b3v07.fa_mod.html -o 1_Tps_b3v07_R_data.txt

```

2) See the R script in R/Repeatmasker dir.

## dnaPipeTE output to R plot (R dir)

1) Results files in:

```
/home/agbast/Data/transfer/dnaPipeTE_landscape_files.zip
```

2) See the R script in R/dnaPipeTE dir.

## polymorphism

similar to previous, just for all reseq
```
map_reads_pol.sh
htseq_poly.sh
count_pol.sh
```



