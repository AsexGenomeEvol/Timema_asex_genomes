### Exploration of SVs

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
source('scripts/vcf_processing_fctions.R')

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

Now, I (need?) want to create a table of filtering thesholds for individual genomes using the figures I just generated, so I make `tables/SVs/coveage_thresholds.tsv` with more relaxed and more stringent thresholds.

Now I want to filter the variants

- no filtering (`diploidSV_reduced.vcf`)
- relaxed: remove only those with split read OR read pair support greater than relaxed threshold (`diploidSV_filt_relaxed.vcf`)
- stringent: split read OR read pair support in the stringent interval (`diploidSV_filt_stringent.vcf`)
- very stringent: split read support in the stringent interval (`diploidSV_filt_very_stringent.vcf`)


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
