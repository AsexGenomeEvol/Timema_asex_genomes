

# F - snp calling (final round) :

(go to previous step : [base recalibration](../E_recalibration))

---------

Base quality scores now being recalibrated, we can perform a definitive round of variant (SNP/indel) calling. Explanations about the pipeline are mainly given here : [snp calling (first round)](../D_snp_calling_round0). In this section, we mostly emphasise on the differences between the two runs.

Following commands are given for species '**1_Tdi**' (samples : *Tdi_01*, *Tdi_02*, *Tdi_03*, *Tdi_04*, *Tdi_05*).



### 1) per-sample calling with 'HaplotypeCaller' :


**HaplotypeCaller** calls SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. 

````
# SAME COMMAND THAN AT FIRST ROUND :
for sample in Tdi_{01,02,03,04,05}
do
java -Xmx20g -jarGenomeAnalysisTK.jar \
   -T HaplotypeCaller \
   -R 1_Tdi_b3v06.fa  \
   -I $sample.bam   \                  # input file
   --genotyping_mode DISCOVERY  \
   --emitRefConfidence GVCF    \
   -o $sample.g.vcf \                  # output file
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
   -V  1_Tdi.gvcf.list  \          # input file
   -allSites \                     # additional option
   -nt 10    \        
   -o  1_Tdi.allVariant_raw.vcf    # output file
````
#### output file :
`1_Tdi.allVariant_raw.vcf`

**note:** we add the `-allSites` tag in order to print **all** positions in the output vcf (including **monomorphic** positions which represent the majority of the positions).
Despite the large file size produced (~25 Go), this is necessary to later calculate some statistics such as heterozygosity (for which we need to apply equivalent criteria to monomorphic sites in order to correctly evaluate the polymorphism); or to convert the vcf files into fasta files (which are then used in some of our pipelines).


### 3) extract SNPs and indels :

Still, we create two additional vcf by extracting only SNP or indel positions (similarly to first round).
````
# this time, we use 'GATK SelectVariants' instead of 'vcftools' for better running time.

# Get SNPs only :
java -jar GenomeAnalysisTK.jar \
    -T SelectVariants          \
    -R 1_Tdi_b3v06.fa          \
    -V 1_Tdi.allSite_raw.vcf   \    # input
    -selectType SNP            \
    -o 1_Tdi.SNP_raw.vcf            # output
# remove deleted positions (that are still written) :
grep -vP "\t\*\t" 1_Tdi.SNP_raw.vcf > 1_Tdi.SNP_raw.vcf.noDel
mv 1_Tdi.SNP_raw.vcf.noDel 1_Tdi.SNP_raw.vcf


# Get indels only :
java -jar GenomeAnalysisTK.jar \
    -T SelectVariants          \
    -R 1_Tdi_b3v06.fa          \
    -V 1_Tdi.allSite_raw.vcf   \    # input
    -selectType INDEL          \
    -o 1_Tdi.indel_raw.vcf          # output
````

#### output files :
````
1_Tdi.SNP_raw.vcf
1_Tdi.indel_raw.vcf
````


### 4) (visually) determine thresholds for hard filtering :

Distributions were very similar between the two rounds, therefore we did not change any threshold value for SNPs or indels, but we added new ones, see below :

#### thresholds for SNPs :
````
# previous thresholds (population level) :
QD < 5.0                       # Quality by Depth (weight call score by the coverage)
FS > 50.0                      # Fisher Strand (probability of strand bias at the site)
SOR > 3.0                      # Strand Odds Ration (another way to estimate strand bias)
MQ < 55.0                      # Mapping Quality 
MQRankSum < -1.0               # do read with ALT allele have a different mapping quality 
ReadPosRankSum < -2.5          # are read positions for REF allele different than ALT allele
# new thresholds (population level) :
DP > 141                       # max total coverage (value differs among species!)
DP < 10                        # min total coverage
# new thresholds (sample level) : 
DP < 10                        # min sample coverage
GQ < 30                        # Genotype Quality (confidence in the sample genotype)
````

#### thresholds for indels :
````

# previous thresholds (population level) :
QD < 5.0
FS > 25.0
SOR > 3.5
ReadPosRankSum < -2.5
# new thresholds (population level) : 
DP < 10                        # min total coverage
# new thresholds (sample level) : 
DP < 10                        # min sample coverage
````

**remarks :**
* **population vs. sample level filters:** in a vcf file, each position has two **FILTER** fields, the first one (which has the `.` value before filtration) concerns the position as a whole, it indicates if the position is variable at the population  scale (its value becomes `PASS` after filtering if it does not fail any test, otherwise it will have a (list of) tag(s) corresponding to failed tests (ex: `lowQD;highFS`)). The second filter field (**FT**) is a per-sample value that indicates if we can trust the genotype of a particular sample (which can be either homozygous or heterozygous) at this position. Therefore, a position can have `PASS` in the **FILTER** field (meaning we are confident there is a SNP at this position), but some specific sample genotypes can still appear dubious (and not getting the `PASS` tag).
* **minimum DP:** variants with low coverage (typically below *10x*) also tend to have low (calling/genotyping) qualities so they have great chances to be filtered out anyway, but to keep things comparable, we also want to filter monomorphic positions that do not have sufficient coverage (as no filter is applied on their quality), hence the use of a minimum coverage filter.
* **maximum DP:** we also added a maximum coverage threshold which is mainly designed to remove positions that could correspond to repeated regions/elements (again, this criterion is also applied to monomorphic sites). The maximum coverage threshold is **species-specific** but always corresponds to : *1.8 x mean(coverage)*, (the mean coverage being calculated on the first 50Mb of the assembly); see below for the coverage distribution in *1_Tdi*. 

![DP](DP_Tdi.png)









#### parameter distributions for SNPs :





### 5) apply hard filters to variants :

A position satisfying any of these criteria will be filtered out.

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

