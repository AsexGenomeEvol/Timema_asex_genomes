
# C - mapping :


Mapping was performed independently for each run (three runs per individual).

Following commands are given for sample **Tdi_01** (run **L2_OBIWAN_314**).


### Pre-mapping :

Assign a read group (**RG**) to the run (it will be added in the bam file headers) :
````
@RG   ID:L2_OBIWAN_312   DT:2017-09-20   LB:lib-L2_OBIWAN_314   PL:ILLUMINA   SM:Tdi_01
````
* **ID**: **read group unique identifier for the run**.
* **LB**: DNA preparation library identifier (used by MarkDuplicates to remove molecular duplicates).
* **PL**: platform/technology used to produce the reads.
* **SM**: **sample** (several read groups can correspond to the same sample).

### Mapping :

#### Command for paired-end fastq files :

````
bwa mem -R $readGroup   \
    -M                  \
    -t 40               \
    $referenceGenome    \
    $fastqR1 $fastqR2   >  $bwaPE.sam
````


#### Command for single-end fastq file :

````
"bwa mem -R $readGroup   \
    -M                   \
    -t 40                \
    $referenceGenome     \
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




