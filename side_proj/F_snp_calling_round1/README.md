

# F - snp calling (final round) :

(go to previous step : [base recalibration](../E_recalibration))

---------

Base quality scores now being recalibrated, we can perform a definitive round of variant (SNP/indel) calling. Explanations about the pipeline are mainly given here : [snp calling (first round)](../D_snp_calling_round0). In this section, we only mostly on the differences between the two runs.

Following commands are given for species '**1_Tdi**' (samples : *Tdi_01*, *Tdi_02*, *Tdi_03*, *Tdi_04*, *Tdi_05*).



### 1) per-sample calling with 'HaplotypeCaller' :


**HaplotypeCaller** calls SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. 

````
for sample in Tdi_{01,02,03,04,05}
do
java -Xmx20g -jarGenomeAnalysisTK.jar \
   -T HaplotypeCaller \
   -R 1_Tdi_b3v06.fa  \
   -I $sample.bam   \
   --genotyping_mode DISCOVERY  \
   --emitRefConfidence GVCF    \
   -o $sample.g.vcf \
   -hets 0.001   \
   -nct  1
 done
````
* `-hets`: prior on heterozygosity level. Its value was set at `0.01` for sexual species and `0.001` for asexuals.


### 2) joint genotyping with 'GenotypeGVCFs' :


'**1_Tdi.gvcf.list' content :**
````
Tdi_01.g.vcf
Tdi_02.g.vcf
Tdi_03.g.vcf
Tdi_04.g.vcf
Tdi_05.g.vcf
````

````
# this step can be parallelised by giving regions to the program
# instead of the complete assembly (as shown below).
java -Xmx20g -jar GenomeAnalysisTK.jar \
   -T  GenotypeGVCFs \
   -R  1_Tdi_b3v06.fa \
   -V  1_Tdi.gvcf.list  \
   -nt 10   \        
   -o  1_Tdi.allVariant_raw.vcf
````


### 3) extract SNPs and indels :


````
# Get SNPs only :
vcftools \
    --vcf 1_Tdi.allVariant_raw.vcf \
    --remove-indels      \
    --out 1_Tdi.SNP_raw  \
    --recode             \
    --keep-INFO-all


# Get indels only :
vcftools \
    --vcf 1_Tdi.allVariant_raw.vcf \
    --keep-only-indels     \
    --out 1_Tdi.indel_raw  \
    --recode               \
    --keep-INFO-all
````

#### output files :
````
1_Tdi.SNP_raw.vcf
1_Tdi.indel_raw.vcf
````

**note:** these two sets of variants are not definitive 
(although some tests were applied by *GATK* to remove the most dubious ones).
The **vcf** files contain the list of possible variant positions (one per line), 
each reported position being associated with a lot of extra-information (such as its global and per-sample coverages, 
and different statistics like the score of the call, mapping quality,...).
It is therefore necessary to set criteria to distinguish between real and false variant (see below).

--------------

The best way to filter the resulting variant callset according to *GATK* recommandations 
is to use *variant quality score recalibration* (**VQSR**),
which uses machine learning to identify annotation profiles of variants that are likely to be real. 
The drawback of this sophisticated method is that it requires a highly curated sets of known variants; 
this makes it impossible to apply it to non-model organisms.
Therefore, in our case, the alternative is to use a quite brutal approach where we will use hard thresholds 
on chosen parameters to manually filter the variants...

#### Some links:
https://gatkforums.broadinstitute.org/gatk/discussion/3225/how-can-i-filter-my-callset-if-i-cannot-use-vqsr-recalibrate-variants
https://software.broadinstitute.org/gatk/documentation/article.php?id=2806



### 4) (visually) determine thresholds for hard filtering :

Here are the parameter threshold values that were chosen for **filtering out** SNPs and indels, 
followed by plots of their distributions (for SNPs) to justify our choices
(although our values differ a bit from those recommended by *GATK* on their own human dataset,
we made use of their graphical recommendations:
https://software.broadinstitute.org/gatk/documentation/article.php?id=6925).

#### thresholds for SNPs :
````
QD < 5.0                       # Quality by Depth (weight call score by the coverage)
FS > 50.0                      # Fisher Strand (probability of strand bias at the site)
SOR > 3.0                      # Strand Odds Ration (another way to estimate strand bias)
MQ < 55.0                      # Mapping Quality 
MQRankSum < -1.0               # do read with ALT allele have a different mapping quality 
ReadPosRankSum < -2.5          # are read positions for REF allele different than ALT allele
````

#### thresholds for indels :
````
QD < 5.0
FS > 25.0
SOR > 3.5
ReadPosRankSum < -2.5
````

#### parameter distributions for SNPs :


![QD](plots/QD_snp.png)

**note:** it is interesting to remark that while sexual species exhibit the expected distribution pattern for the **QD** parameter, asexual species (in *grey*, *black*, *blue*, *green* and *violet*) seem to have a different shape, mostly lacking the central pic (the pic on the left corresponding to errors).

![FS](plots/FS_snp.png)
![SOR](plots/SOR_snp.png)
![MQ](plots/MQ_snp.png)
![MQRankSum](plots/MQRankSum_snp.png)
![ReadPosRankSum](plots/ReadPosRankSum_snp.png)



### 5) apply hard filters to variants :

#### command for SNP vcf :
````
java -jar GenomeAnalysisTK.jar \
    -T VariantFiltration    \
    -R 1_Tdi_b3v06.fa       \
    -V 1_Tdi.SNP_raw.vcf    \
    -o 1_Tdi.SNP_filter.vcf \
    --filterExpression  "QD < 5.0"    --filterName "badQD"  \
    --filterExpression  "FS > 50.0"   --filterName "badFS"  \
    --filterExpression  "SOR > 3.0"   --filterName "badSOR" \
    --filterExpression  "MQ < 55.0"   --filterName "badMQ"  \
    --filterExpression  "MQRankSum < -1.0"          --filterName "badMQRankSum"    \
    --filterExpression  "ReadPosRankSum < -2.5"     --filterName "badReadPosRankSum"  
````

#### output :
````
1_Tdi.SNP_filter.vcf
1_Tdi.indel_filter.vcf
````
These first sets of variants will be used in the following **BQSR** step to mask their position (so they are not mistaken for sequencing errors).

