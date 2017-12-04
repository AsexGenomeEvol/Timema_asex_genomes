
# B - read trimming :


Paired-end reads were trimmed for presence of adapters and quality using **Trimmomatic**.
Similar parameters to the individual used for genome assembly were employed ([Illumina adapters](./Illumina-PE_adapters.fa) used).

Command example (for a single run of a single sample, R1&R2 files processed together) :

````
java -jar trimmomatic-0.36.jar PE -threads 25 \
Tdi_01_L2_OBIWAN_312_R1.fastq.gz         Tdi_01_L2_OBIWAN_312_R2.fastq.gz            \
Tdi_01_L2_OBIWAN_312_R1.cleaned.fastq.gz Tdi_01_L2_OBIWAN_312_R1.se.cleaned.fastq.gz \
Tdi_01_L2_OBIWAN_312_R2.cleaned.fastq.gz Tdi_01_L2_OBIWAN_312_R2.se.cleaned.fastq.gz \
ILLUMINACLIP:Illumina-PE_adapters.fa:3:25:6 LEADING:9 TRAILING:9 SLIDINGWINDOW:4:15 MINLEN:96
````

#### Output files :

Surviving pairs of reads :
````
Tdi_01_L2_OBIWAN_312_R1.cleaned.fastq.gz`
Tdi_01_L2_OBIWAN_312_R2.cleaned.fastq.gz
````
The two output fastq files containing unpaired reads :
````
Tdi_01_L2_OBIWAN_312_R1.se.cleaned.fastq.gz`
Tdi_01_L2_OBIWAN_312_R2.se.cleaned.fastq.gz
````
were merged together into :
````
Tdi_01_L2_OBIWAN_312.se.cleaned.fastq.gz`
````


#### Results:

Around 10% of the reads are trimmed for each run, see detailed results [here](./trimmomatic.log) and final number of reads 
per file [here](./number_reads.csv).
