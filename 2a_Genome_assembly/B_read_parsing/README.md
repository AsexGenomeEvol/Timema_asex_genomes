## Testing and details

We tested two tools for delinking (`Skewer` and `NxTrim`) by comparing distributions of insert sizes of delinked reads mapped back to a draft assembly.
`NxTrim` was selected since it managed to produce less errors in linker identification (evaluated by visual comparison of distributions produced by script `n2_insert_sizes_distributions.R`, insert sizes were computed from mapping files using scipt `n1_get_insert_sizes.py`). I found out that the best will be just to use all mate-pairs and unknown reads from NxTrim.

The correction of reads using `bcf`, but correction of reads was latter shown to decrease quality of assembly, therefore this step was skipped. The script used for correction is `n3_bcf_lsf.sh`.


#### Reevaluation of mate pair reads

following section describes testing, it does not describe any result which would be directly used in the flow. It just justifies some software / parameter choices.
Scripts used for testing are prefixed by `n` (applies to whole repository).

lot of mate-apirs are seqeunced without found linking adapter, those sequences can be both contaminants or true mate-pairs. I will try software Skewer and Nxtrim to trim a filter mate-pairs. Skewer seems to be more relaxing in accepting mate pairs and also to be better in adapter trimming. NxTrim is more stringent in filtering (produces less mate-pairs from the same input), doing worse job in trimming, but it is able to categorise reads to 4 categories - se, pe, mp and unknown. I will map back pe, mp and unknown reads of nxtrim and mp of skewer to find how many conflict to our assembly I will find.

The tests will be run on Tdi and Tps (Darren's favourire).

##### Trimming step

NxTrim

```sh
module add UHTS/Assembler/NxTrim/0.4.1
nxtrim -1 160203_SND405_B_L004_HYI-31_R1.fastq.gz -2 160203_SND405_B_L004_HYI-31_R2.fastq.gz -O is_3000 --preserve-mp --separate
```

Skewer

```sh
/home/kjaron/bin/skewer -e -L 120 -n -q 26 -l 24 -t 8 -m mp -o Tdi_t2 160203_SND405_B_L004_HYI-31_R1.fastq.gz 160203_SND405_B_L004_HYI-31_R2.fastq.gz
```

##### Mapping

was performed using bwa-mem. First of all, the index of a reference has to be build.
All already computed indexes will be stored with assemblies in `/scratch/beegfs/monthly/kjaron/timema_asm/*`

```
bwa index <genome.fa>
```

Finally the script submiting a job to Vital-it has following syntax

```
bash scripts/map_pair_end_lsf.sh <output> <reference_with_index> <R1.fq> <R2.fq>
```

mapping of skewer and nxtrim reads was done using `scripts/map_pair_end_lsf.sh`.
All the output is written in `.bam` files and it was parsed using `pysam`.

##### Insert sizes

will be computed using a python scipt `n1_get_insert_sizes.py`. it is so far working on `.sam` file sorted by name. Insert sizes are computed from primary mappings only. Tab (space, to be fixed) delimited table is printed on the standard output with 4 columns: `pair_id mp mp_is pe_is`, where `mp` is 1 is the orientation is mate-pair ( <-- --> ), 0 if the orientation is pe (  --> <-- ) and -1 if orientations is unknown, mp_is and `pe_is` are minimal mate-pair and pair-end distance if contigs where reads map would be "glued" together by one or the other side. `Rscript n2_insert_sizes_distributions.R` is used to plot log histograms of insert sizes, two arguments are passed, table of insert sizes and a name of the histogram to plot.

```sh
samtools sort -n -O 'sam' Tdi_nxt_3k_unknown.bam Tdi_nxt_3k_unknown.sorted.sam
python n1_get_insert_sizes.py Tdi_nxt_3k_unknown.sorted.sam > Tdi_is_unknown_reads.tsv
Rscript n2_insert_sizes_distributions.R Tdi_is_unknown_reads.tsv Tdi_is_unknown.pdf
```
Analogically, the rest of log histograms can be ploted (files `*.pdf`).

##### Conclusion

Conclusion here is, that `NXtrim` makes a good job how it is and remapping reads on draft assembly is not making it any better. Also seems that "unknown" category is just as good as "mp" category. Overlapping reads have loads of very very short reads (therefore those pairs are not overlapping at all), I have for allpaths-lg filtered too short reads by

```
sickle pe -f ../Tdi_R1t_is_225.fq -r ../Tdi_R2t_is_225.fq -t sanger -o Tdi_R1tf_is_225.fq -p Tdi_R2tf_is_225.fq -s Tdi_se.fq -l 80
```

Skewer is apparently confused and it contains a lot of pair-end reads (1 / 4).
NXtrim does a better job, mate pairs have expected insert size of reads mapping
to one scaffold only.  However 1 / 10 of them in wrong orientation. However,
there is not apriory way how to distinguish them. I would have to assume that the
assembly is correct or that sequencing gone wrong or that both NXtrim and skewer
are messing with orientation (unlikely). Category `unknown` have completely same
distribution as mate-pairs. Therefore I will use them for scaffolding, but not
for error correction afterwards. pair-end reads of NXtrim are pretty robust.

Long story short, trust NXtrim!

#### Short read correction

I will try to correct reads of `Tcm` on dee-server04 to find resources needed for read
correction and at the same time to find out if the assembly is going to improve with
the read correction.

 ```
mkdir -p /scratch/local/kjaron/temp /scratch/local/kjaron/2_Tcm_corrected_reads
export TMP_DIR=/scratch/local/kjaron/temp
cd /scratch/local/kjaron/2_Tcm_corrected_reads
/usr/bin/time -f '%M %E %P' bfc -s 1.3g -t16 $TROOT/data/2_Tcm/trimmed_reads/Tcm_R1t_is_225.fq.gz | gzip -1 > Tcm_R1tc_is_225.fq.gz
 ```

for this dataset it took 20G of memory in peak. CPU: 1729.746 sec is stable (i.e. multiuthreading is effective).
To estimate requirements I will run the biggest file I have as well - CPU: 26766.052 sec; mem peak: 40G. CPU time goes linear
with file size: 1.5 CPU minutes per 1G of file.

According some tests (that should be described in secton C), corrected reads assemble better than uncorrected reads, therefore I correct all reads I got:

```bash
cd $TROOT/data
for SP in $(echo 1_T*); do
	mkdir -p $SP/corrected_reads
	cd $SP/trimmed_reads
	READ_FILES=$(ls T*)
	# I want logs to be saved in this folder
	cd ../corrected_reads
	for READ_FILE in $READ_FILES; do
		$TROOT/B_read_parsing/6_bcf_lsf.sh $SP $READ_FILE
	done
	cd $TROOT/data
done
```
According some more tests, it actually makes things worse, previous positive results were caused by adding more data at the same time as introduction of corrected reads, therefore I will NOT use corrected reads by bcf and I renamed script to `n3_bcf_lsf.sh`.

### kmergenie

We need to predict a kmer for the genome assembly.
According manual of kmergenie, a tool for predicition of optimal assembly kmer, a parameter for the diploid mode of program should be specified only in cases that heterozygosity of the organism is high. Since I knew, that heterozygosity is not in any way extremely high I simply run the program with the default setting with restriction of kmer choices given our coverage :

```
kmergenie read_files.list -l 59 -k 89 -s 6 -t 16
```

where `read_files.list` is a list of all sequencing reads on the input of genome assembly.

HOWEVER, Patrick have run the same program with `--diploid` flag and his kmers have lead to better assemblies than mine. We will not reassemble all the genomes again, however we should keep this in mind if for any reason we decide to do so.

### practical notes

- Script for trimming with trimmomatic was removing non-empty directory, therefore all logs are with exit code 1, even though they are fine.
