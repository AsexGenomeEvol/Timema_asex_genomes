### SV calling

Complementing nucleotide-based measures we searched for **Structural variants**, to investigate if any heterozygosity is carried by SVs. We used Manta (v1.5.0), a diploid-aware pipeline for structural variant (SV) calling, in the same set of re-sequenced individuals used for SNP heterozygosity estimates.

We found a high frequency of heterozygous SVs with approximately twice the expected coverage (SM Figure 7), which likely represent false positives. To reduce the number of false positives, we filtered very short SVs (30 bases or less) and  kept only variant calls that had either split read or paired-end read support within the expected coverage range, where the coverage range was defined individually for each sample by manual inspection of coverage distributions. The filtered SV calls were subsequently merged into population SV calls using SURVIVOR (v1.0.2). The merging criteria were: SV calls of the same type on the same strand with breakpoints distances shorter than 100 bp.  



 -  [scripts](scripts)

 -  [C_variant_analysis](C_variant_analysis)


Alright, all the SV calls take for ever detecting breakpoint between scaffolds as SVs (scaffolds are considered chromosomes in the SV world). To prevent this unintended behaviour I need to get rid of all reads mapping to edges of different scaffolds.
I was considering `BAMQL`, but it was too complicated to install on the cluster.
So, I will map `A_mapping_and_preprocessing/filter_splitreads.py`

## Structural variations

Detection of structural variants from population data will be based on pair-end illumina data. We are interested in difference between sexual and asexual species, we therefore use a SV caller well balancing specificity and sensitivity [Manta](https://github.com/Illumina/manta), due to the problem of collapsed duplicates we also introduce coverage, based filtering and use a [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR) to merge calls (variants in different individuals of the same type with breakpoints less than 100 bases from each other).

#### SV calls: pair-end

Indexing mapped reads to the reference gneomes (`data/final_references/*_b3v08.fasta.gz`)

```
for bam in data/mapped_reseq_reads/*.bam; do
  qsub -o logs/ -e logs/ -cwd -N bamindex -V -pe smp64 1 -b yes "samtools index $bam"
done
```

I downloaded precompiled binaries of (Manta)[https://github.com/Illumina/manta] v1.0.3. Using commands in the script `E_structural_variation_calling/manta.sh` I get 50 manta SV calls: `data/manta_SV_calls/data/"$SP"/"$IND"_manta/results/variants/diploidSV.vcf`

#### SV qc and filtering before merging

Removing header from vcf files (non-necessary for filtering)

```
for file in $(find ./data/manta_SV_calls -name "diploidSV.vcf.gz"); do
    out=${file%.vcf.gz}_reduced.vcf
    echo $file $out
    zcat $file | grep -v "^#" > $out
done
```

in the following R snippet I will plot all the coverage distributions of split reads and pair end coveages of homozygous and heterozygous variants per individual.

```{R}
library(AsexStats)
source('E_structural_variation_calling/vcf_processing_fctions.R')

sp <- timemas$codes[1]

for(sp in timemas$codes){
    sp_short <- substr(sp, 3, 5)

    SV_files <- paste0('data/manta_SV_calls/data/', sp, '/', sp_short, '_0', c(1:5) ,'_manta/results/variants/diploidSV_reduced.vcf')

    pdf(paste0('figures/manta_SV_coverages_', sp ,'.pdf'))
        for (ind in 1:5){
            SV_tab <- read.table(SV_files[ind], header = F)
            print(plot_one(SV_tab, 'SR', paste0(sp_short, '0', ind, ' SR')))
            print(plot_one(SV_tab, 'PR', paste0(sp_short, '0', ind, ' PR')))
        }
    dev.off()
}
```

Now, I (need?) want to create a table of filtering thesholds for individual genomes using the figures I just generated, so I make `tables/SVs/coveage_thresholds.tsv` with more relaxed and more stringent thresholds. The following snippet will generate a few filtering files based on different stringency levels.

```R
library(AsexStats)
source('E_structural_variation_calling/vcf_processing_fctions.R')

filtering_thresholds <- read.table('tables/SVs/coveage_thresholds.tsv', header = T, row.names = 1)

for(sp in timemas$codes){
    sp_short <- substr(sp, 3, 5)

    SV_files <- paste0('data/manta_SV_calls/data/', sp, '/', sp_short, '_0', c(1:5) ,'_manta/results/variants/diploidSV_reduced.vcf')
    out_relaxed <- paste0('data/manta_SV_calls/data/', sp, '/', sp_short, '_0', c(1:5) ,'_manta/results/variants/diploidSV_filt_relaxed.vcf')
    out_str <- paste0('data/manta_SV_calls/data/', sp, '/', sp_short, '_0', c(1:5) ,'_manta/results/variants/diploidSV_filt_stringent.vcf')
    out_very_str <- paste0('data/manta_SV_calls/data/', sp, '/', sp_short, '_0', c(1:5) ,'_manta/results/variants/diploidSV_filt_very_stringent.vcf')

    for (ind in 1:5){
        SV_tab <- read.table(SV_files[ind], header = F)
        SR_cov_1 <- sapply(strsplit(SV_tab[,10], ':'), str2depth, 6, 1)
        SR_cov_2 <- sapply(strsplit(SV_tab[,10], ':'), str2depth, 6, 2)
        SR_cov <- SR_cov_1 + SR_cov_2
        PR_cov_1 <- sapply(strsplit(SV_tab[,10], ':'), str2depth, 5, 1)
        PR_cov_2 <- sapply(strsplit(SV_tab[,10], ':'), str2depth, 5, 2)
        PR_cov <- PR_cov_1 + PR_cov_2

        # relaxed: remove only those with split read OR read pair support greater than relaxed threshold
        relaxed_u <- filtering_thresholds[paste0(sp_short,'0',ind), 'relaxed_u']
        keep <- is.na(SR_cov) | SR_cov < relaxed_u
        write.table(SV_tab[keep,], out_relaxed[ind], quote = F, sep = '\t', row.names = F, col.names = F)

        stringent_l <- filtering_thresholds[paste0(sp_short,'0',ind), 'stringent_l']
        stringent_u <- filtering_thresholds[paste0(sp_short,'0',ind), 'stringent_u']

        # stringent: split read OR read pair support in the stringent interval
        SR_in_range <- SR_cov > stringent_l & SR_cov < stringent_u
        SR_in_range[is.na(SR_in_range)] <- F # missing value means 0
        PR_in_range <- (PR_cov > stringent_l & PR_cov < stringent_u)
        PR_in_range[is.na(PR_in_range)] <- F # missing value means 0
        write.table(SV_tab[SR_in_range | PR_in_range,], out_str[ind], quote = F, sep = '\t', row.names = F, col.names = F)

        # very stringent: split read support in the stringent interval
        write.table(SV_tab[SR_in_range,], out_very_str[ind], quote = F, sep = '\t', row.names = F, col.names = F)
    }
}
```

Now we have variants:

- no filtering (`diploidSV_reduced.vcf`)
- relaxed: remove only those with split read OR read pair support greater than relaxed threshold (`diploidSV_filt_relaxed.vcf`)
- stringent: split read OR read pair support in the stringent interval (`diploidSV_filt_stringent.vcf`)
- very stringent: split read support in the stringent interval (`diploidSV_filt_very_stringent.vcf`)


##### SV curation

This is sort of a post-analysis step to make sure what we do makes sense.

We will use [SV-plaudit](https://github.com/jbelyeu/SV-plaudit), a high throughput SV curation tool. It internally uses [samplot](https://github.com/ryanlayer/samplot) to plot individual SVs, I will try to use it first internally, then I will run the whole framework

```
samplot plot -n Tms_00 Tms_01 Tms_02 Tms_03 Tms_04 Tms_05 -b data/mapped_reseq_reads/Tms_00_to_b3v08_mapped_within_scfs.bam data/mapped_reseq_reads/Tms_01_to_b3v08_mapped_within_scfs.bam data/mapped_reseq_reads/Tms_02_to_b3v08_mapped_within_scfs.bam data/mapped_reseq_reads/Tms_03_to_b3v08_mapped_within_scfs.bam data/mapped_reseq_reads/Tms_04_to_b3v08_mapped_within_scfs.bam data/mapped_reseq_reads/Tms_05_to_b3v08_mapped_within_scfs.bam -o data/sandbox/3_Tms_b3v08_scaf000092_INV_80652_80811.png -s 80652 -e 80811 -c 3_Tms_b3v08_scaf000092 -a -t INV
```

this looks bad. Like really bad. Let's try to get them for all

```
qsub -o logs/ -e logs/ -cwd -N samplot -V -pe smp64 1 -b yes 'samplot vcf --vcf data/genotyping/3_Tms_merged_calls_naive_header.vcf -d figures/SVs -O png -b data/mapped_reseq_reads/Tms_00_to_b3v08_mapped_within_scfs.bam data/mapped_reseq_reads/Tms_01_to_b3v08_mapped_within_scfs.bam data/mapped_reseq_reads/Tms_02_to_b3v08_mapped_within_scfs.bam data/mapped_reseq_reads/Tms_03_to_b3v08_mapped_within_scfs.bam data/mapped_reseq_reads/Tms_04_to_b3v08_mapped_within_scfs.bam data/mapped_reseq_reads/Tms_05_to_b3v08_mapped_within_scfs.bam --sample_ids Tms_00 Tms_01 Tms_02 Tms_03 Tms_04 Tms_05 > F_structural_variation_analysis/samplot_commands.sh'
```

It seems that all heterozygous SVs in asexuals are shaky.

##### SVs in a Long Read asexual individual

We have one `Tdi` individual, it's long reads are `data/1_Tdi/pacbio_1Tdi/raw_reads/`. We will use a classical nglmr/sniffles SV pipeline.

```
zcat data/1_Tdi/pacbio_1Tdi/raw_reads/*.fastq.gz | ngmlr -t 16 -r data/final_references/1_Tdi_b3v08.fasta | samtools view -bh - > data/1_Tdi/pacbio_1Tdi/mapped_to_b3v08.bam
samtools index data/1_Tdi/pacbio_1Tdi/mapped_to_b3v08.bam
sniffles -m data/1_Tdi/pacbio_1Tdi/mapped_to_b3v08.bam -v data/1_Tdi/pacbio_1Tdi/Tdi06_sniffles.vcf
```

This took for ever! Either I must filter out the bamfile out of all "between scaffolds" mapping, but Tanja suggested that I could just make quick and dirty long read assembly and call SVs on that.

Quick and dirty redbeans assembly

```
conda activate /ceph/users/chodson/.conda/envs/pbassembly

qsub -o logs/ -e logs/ -cwd -N read_len_dist -V -pe smp64 1 -b yes '~/generic_genomics/fastq2read_lengths.sh data/1_Tdi/pacbio_1Tdi/raw_reads/*.fastq.gz > data/1_Tdi/pacbio_1Tdi_readlen_dist.tsv'

qsub -o logs/ -e logs/ -cwd -N redbeans -V -pe smp64 32 -b yes 'wtdbg2 -L 1000 -x preset3 -t 32 -g 1300m -i data/1_Tdi/pacbio_1Tdi/raw_reads/7k1yOJfO_m54161_171122_021109.fastq.gz -i data/1_Tdi/pacbio_1Tdi/raw_reads/FHbycsRw_m54161_190117_061402.fastq.gz -i data/1_Tdi/pacbio_1Tdi/raw_reads/tW3dEm3W_m54161_190213_054535.fastq.gz -i data/1_Tdi/pacbio_1Tdi/raw_reads/u2gBb1JT_m54161_190215_191639.fastq.gz -i data/1_Tdi/pacbio_1Tdi/raw_reads/Y1EnuWX5_m54161_171127_212638.fastq.gz -o data/1_Tdi/pacbio_1Tdi/redbeans_asm/Tdi_pb_asm00 && wtdbg-cns -t 32 -i data/1_Tdi/pacbio_1Tdi/redbeans_asm/Tdi_pb_asm00.ctg.lay.gz -o data/1_Tdi/pacbio_1Tdi/redbeans_asm/Tdi_pb_asm00.ctg.fa'
```

and map the reads on the new reference.

```
qsub -o logs/ -e logs/ -cwd -N sniffles -V -pe smp64 32 -b yes "zcat data/1_Tdi/pacbio_1Tdi/raw_reads/*.fastq.gz | ngmlr -t 32 --no-progress -r data/1_Tdi/pacbio_1Tdi/redbeans_asm/Tdi_pb_asm00.ctg.fa | samtools view -bh - > data/1_Tdi/pacbio_1Tdi/mapped_to_pb_asm00.bam"
```

run sniffles

```
qsub -o logs/ -e logs/ -cwd -N sniffles -V -pe smp64 32 -b yes "mkdir -p /scratch/kjaron/snif_temp && sniffles -t 32 -q 15 -m data/1_Tdi/pacbio_1Tdi/mapped_to_pb_asm00.bam -v data/1_Tdi/pacbio_1Tdi/Tdi06_on_pb_asm_sniffles.vcf --tmp_file /scratch/kjaron/snif_temp && rm -r /scratch/kjaron/snif_temp"
```

and visualise them:

```
qsub -o logs/ -e logs/ -cwd -N sortbam -V -pe smp64 16 -b yes 'mkdir -p /scratch/kjaron && samtools sort -T /scratch/kjaron/sorting_bam -o data/1_Tdi/pacbio_1Tdi/mapped_to_pb_asm00_sorted.bam -@16 data/1_Tdi/pacbio_1Tdi/mapped_to_pb_asm00.bam && samtools index data/1_Tdi/pacbio_1Tdi/mapped_to_pb_asm00_sorted.bam'
```

```
qsub -o logs/ -e logs/ -cwd -N samplot -V -pe smp64 1 -b yes 'START=3019763; END=3019890; SCF=ctg152; TYPE=DEL; samplot plot -n Tdi_PB -b data/1_Tdi/pacbio_1Tdi/mapped_to_pb_asm00_sorted.bam -o data/sandbox/"$SCF"_"$START"_"$END"_"$TYPE".png -s $START -e $END -c $SCF -a -t $TYPE'
```



```
while read SCF START END TYPE
do
    echo samplot plot -n Tdi_PB -b data/1_Tdi/pacbio_1Tdi/mapped_to_pb_asm00_sorted.bam -o data/sandbox/"$SCF"_"$START"_"$END"_"$TYPE".png -s $START -e $END -c $SCF -a -t $TYPE >> plotting_commands.sh
done < tables/Tdi_PB_SVs.tsv
```

```
parallel -j 1 'qsub -o logs/ -e logs/ -cwd -N samplot -V -pe smp64 1 -b yes {}' :::: plotting_commands.sh
```

```
qsub -o logs/ -e logs/ -cwd -N samplot -V -pe smp64 1 -b yes 'samplot vcf --vcf data/1_Tdi/pacbio_1Tdi/Tdi06_on_pb_asm_sniffles.vcf -d figures/SVs -O png -b data/1_Tdi/pacbio_1Tdi/mapped_to_pb_asm00.bam --sample_ids Tdi_PB > F_structural_variation_analysis/samplot_commands_PB.sh'
```

**SVs**

In `data/manta_SV_calls/data/$SP/` directories we have merged files: split read OR read pair support in the stringent interval (`SVs_filt_stringent_union.vcf`)


### Plotting heteoryzogisty

```
D_variant_analysis/make_heterozygosity_table.R
D_variant_analysis/plot_heterozygosity.R
```


TODO

```
5_Heterozygosity/scripts/coverages_table.R
5_Heterozygosity/scripts/prepare_coverage_tables.R
```
