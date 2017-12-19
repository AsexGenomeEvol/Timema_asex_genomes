
# E - base quality score recalibration (BQSR) :

(go to next step : [snp calling (last round)](../F_snp_calling_round1), previous step : [snp calling (first round)](../D_snp_calling_round0))

---------

This step aims to detect systematic errors made by the sequencer when it estimates the quality score of each base call,
and then readjust their values in the **bam** files. This process is accomplished by analyzing the covariation among several features of a base, mainly its *reported quality score*, its *position* within the read and the effect of preceding and current nucleotide (*sequencing chemistry effect*).

#### Some links :
https://software.broadinstitute.org/gatk/documentation/article.php?id=44
https://gatkforums.broadinstitute.org/gatk/discussion/2801/howto-recalibrate-base-quality-scores-run-bqsr

As for the mapping step, the procedure is done independently for each sample. 
Following commands are given for sample : '**Tdi_01**' (species: **1_Tdi**).


### 1) analyze patterns of covariation in the sequence dataset :

````
java -Xmx20g -jar GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R 1_Tdi_b3v06.fa   \
    -I Tdi_01.bam       \
    -knownSites 1_Tdi.SNP_filter.vcf   \
    -knownSites 1_Tdi.indel_filter.vcf \
    -nct 20             \
    -o Tdi_01.recal.table
````
#### output :
`Tdi_01.recal.table` : it contains the covariation data that will be used in the step to recalibrate base qualities.



### 2) do a second pass to analyze covariation remaining after recalibration :


````
java -Xmx20g -jar GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R 1_Tdi_b3v06.fa   \
    -I Tdi_01.bam       \
    -knownSites 1_Tdi.SNP_filter.vcf   \
    -knownSites 1_Tdi.indel_filter.vcf \
    -BQSR Tdi_01.recal.table          \
    -nct 20             \
    -o Tdi_01.post_recal.table
````
#### output :
`Tdi_01.post_recal.table` : second report which will be used in next step to evaluate the effect of base quality recalibration.



### 3) generate before/after recalibration plots :


````
java -Xmx10g -jar GenomeAnalysisTK.jar \
    -T AnalyzeCovariates               \
    -R 1_Tdi_b3v06.fa                  \
    -before Tdi_01.recal.table         \
    -after  Tdi_01.post_recal.table    \
    -plots  Tdi_01.recal_plots.pdf
````
#### output :
`Tdi_01.recal_plots.pdf` : plots that show how the reported base qualities match up to the empirical qualities 
calculated by **BaseRecalibrator** (download [pdf file for Tdi_01 sample](Tbi_01.recal_plots.pdf)). Comparing the before and after plots allows to check the effect of the base recalibration process before applying it to the sequence data (*bam* file).



### 4) apply the recalibration to sequence data :

````
java -jar GenomeAnalysisTK.jar \
    -T PrintReads              \
    -R 1_Tdi_b3v06.fa          \
    -I Tdi_01.bam              \
    -BQSR Tdi_01.recal.table   \
    -nct 20                    \
    -o Tdi_01.recal.bam
````
#### output :
`Tdi_01.recal.bam` : recalibrated *bam* file which will be used to perform the [final round of variant calling](../F_snp_calling_round1).



