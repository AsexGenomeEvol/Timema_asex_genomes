
# C - mapping :


Mapping was performed independently for each run (three runs per individual).

Following commands are given for sample **Tdi_01** (run: **L2_OBIWAN_314**).


### 0) preliminary step :

Assign a read group (**RG**) to the run (it will be added in the bam file headers) :
````
@RG   ID:L2_OBIWAN_312   DT:2017-09-20   LB:lib-L2_OBIWAN_314   PL:ILLUMINA   SM:Tdi_01
````
* **ID**: **read group unique identifier for the run**.
* **LB**: DNA preparation library identifier (used by MarkDuplicates to remove molecular duplicates).
* **PL**: platform/technology used to produce the reads.
* **SM**: **sample** (several read groups can correspond to the same sample).

### 1) mapping :

#### Command for paired-end fastq files :

````
bwa mem -R $readGroup   \
    -M                  \
    -t 40               \
    $referenceGenome_b3v04 \
    $fastqR1 $fastqR2   >  $bwaPE.sam
````


#### Command for single-end fastq file :

````
bwa mem -R $readGroup   \
    -M                   \
    -t 40                \
    $referenceGenome_b3v04  \
    $fastqUnpaired       >  $bwaSE.sam
````


#### Default parameter values :

````
#  -k INT     minimum seed length [19]
#  -w INT     band width for banded alignment [100]
#  -d INT     off-diagonal X-dropoff [100]
#  -r FLOAT   look for internal seeds inside a seed longer than {-k} * FLOAT [1.5]
#  -c INT     skip seeds with more than INT occurrences [500]
#  -A INT     score for a sequence match [1]
#  -B INT     penalty for a mismatch [4]
#  -O INT     gap open penalty [6]
#  -E INT     gap extension penalty; a gap of size k cost {-O} + {-E}*k [1]
#  -L INT     penalty for clipping [5]
#  -U INT     penalty for an unpaired read pair [17]
#  -T INT     minimum score to output [30]
#  -M         mark shorter split hits as secondary (for Picard/GATK compatibility)
#             manual: 'the BWA-MEM algorithm performs local alignment. 
                       It may produce multiple primary alignments for different part of a query sequence. 
                       This is a crucial feature for long sequences. 
                       However, some tools such as Picard markDuplicates does not work with split alignments. 
                       One may consider to use option -M to flag shorter split hits as secondary.'
````

**NOTE:** mapping was not performed on 'final' assemblies (version: **b3v06**), but on a previous version (**b3v04**) containing **all** scaffolds (even those of size <1kb that were then removed in last version). This is to avoid that resequenced reads map on an incorrect scaffold because of the absence of their true target. A quick comparison between mapping results of the same run on both assemblies shows that the amount of concerned reads can be large (up to 20% of the total pool of reads can map incorrectly to the wrong scaffold when then they in fact belong elsewhere). Keep in mind however that such hits are likely to receive a poor **mapping quality score** that will later reduce their influence in the discrimination between true and false SNPs.


### 2) filter sam files on regions and mapping quality :

Only reads mapped to scaffolds present in the final assembly (version: **b3v06**) and above a certain mapping quality threshold (phred-score: **20**) were kept (this also removes unmapped reads) :
````
samtools view -h -q 20 -L $bed -@ 10 $bwaPE.sam > $bwaPE.filtered.sam
samtools view -h -q 20 -L $bed -@ 10 $bwaSE.sam > $bwaSE.filtered.sam
````
`$bed`: bed file containing list of coordinates of final scaffolds.


### 3) convert sam file to sorted bam :

````
samtools sort -@ 10 -O bam -l 9 -T $bwaPE.temp.bam -o $bwaPE.bam $bwaPE.filtered.sam
samtools sort -@ 10 -O bam -l 9 -T $bwaSE.temp.bam -o $bwaSE.bam $bwaSE.filtered.sam
````

### 4) merging paired-end and single-end bam files :

````
samtools merge -f -@ 10 -c $bwa.bam $bwaPE.bam $bwaSE.bam      
````
**note:** `-c` tag is use to keep the same '**@RG:ID**' when it is the same in the different input files (by default, the program adds a suffix to one of the read group!).


### 5) rewrite header section by keeping only scaffolds from final assembly (and indexing) :

Later programs will check the correspondence between the list of scaffolds in the assembly (fasta) and the headers of the bam file (where all scaffolds of **b3v04** assembly are currently present).
````
java -jar picard.jar ReorderSam \
    I=$bwa.bam                  \
    O=$bwa.reordered.bam        \
    R=$referenceGenome_b3v06    \
    S=true                      \
    CREATE_INDEX=TRUE
# replace old bam by new :
mv $bwa.reordered.bam $bwa.bam
````

### 6) mark duplicates with Picard tools (and reindex) :

Duplicates (ie, reads having exactly the same sequence and therefore mapping to the same location), mostly arising from PCR amplification bias introduced during library construction (or corresponding to optical duplicates), are removed from the bam file, as they inflate the coverage of a position and can lead for instance to errors looking like true snps; but more generally, they will lead to a breakdown of the statistical models for variant calling that assume some sort of independence between measurements.

````
java -Xmx25g -jar picard.jar MarkDuplicates \
   INPUT=$bwa.bam \
   OUTPUT=${bwa}_md.bam \
   METRICS_FILE=log/${bwa}_duplicate_metrics.txt \
   MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
   REMOVE_DUPLICATES=true               
mv ${bwa}_md.bam $bwa.bam
# reindexing :
samtools index $bwa.bam
````


### 7) merge all (three) bam files from a single sample :

(after previous steps are done for each of the three runs of a sample)

````
samtools merge -@ 10 -b $sample.run_list $sample.bam        # merge bam files (list in $sample.run_list, output: $sample.bam)
samtools index $sample.bam                                  # reindex
````

**BILAN:** We get a single bam file per sample (5 individuals *x* 10 species = 50 bam files).
