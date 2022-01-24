### Localisation of color loci from literature

Regions in _Timema_ genomes associated with color loci:

_T. cristinae_:

> The major locus associated with color is called _Mel-Stripe_ and it is on the LG8. It spans about 10.5Mbp. The recombination is suppressed possibly due to large scale inversion (Nosil et al. 2018; Doi: 10/gc28jr) and one edge of the locus exhibits a large-scale (~1 megabase pair) deletion (in the green morph, Villoutreix et al. 2020).

Relative coordinates from Figure 1 in Nosil et al. 2018:
 - within 702.1 scaffold: 0.521 - 1
 - within 128 scaffold: 0 - 0.423

_T. californicum_:
The relative position of the deletion within the color locus (the one above): 0.18 - 0.58
Absolute position on the figure: 5.16 - 5.69 Mbp
From legend of Fig. 2:

 > (indel, 5.0 to 6.0 megabase pairs on scaffold 128; Mel-Stripe excluding indel, 1 to 4,139,489 base pairs on scaffold 702 and 1 to 5.0 megabase pairs on scaffold 128).

I need to recalculate it to nucleotides using scaffold lenghts and look at the position of those nt’s in the gneome (I just had the figure in front of me, so I measured it relativelly)

#### calculating values from figures

```{R}
reference <- read.table('data/external_ref/sex_lg_assigment_scores_1.4a.tsv')
colnames(reference) <- c('scf_o', 'scf', 'score', 'cov', 'len', 'asignment')
reference$chromosome <- sapply(strsplit(reference$scf, "_"), function(x) { x[1] } )
lg8 <- reference[reference$chromosome == 'lg8', ]

row.names(reference) <- reference$scf

round(reference['lg8_ord14_scaf702.1', 'len'] * c(0.521, 1))

round(reference['lg8_ord15_scaf128', 'len'] * c(0, 0.423))

```



### Mapping Tcm deletion

```r
library(AsexStats)
library(RColorBrewer)

mapping_file <- paste0('data/b3v08_anchoring_to_LGs/2_Tcm_scf_block_alignment.tsv')
mapping_table <- read.table(mapping_file, header = T)
nrow(mapping_table)

# filter out blocks mapped "unequally" - reference block / query block can not be < 0.5
mapping_table <- mapping_table[abs(mapping_table$block_r_end - mapping_table$block_r_start) / mapping_table$block_size > 0.5, ]
nrow(mapping_table)

# keep only lg8 blocks
mapping_table <- mapping_table[mapping_table$lg == 'lg8', ]
nrow(mapping_table)

mapping_table <- mapping_table[order(mapping_table$lg, mapping_table$lg_start), ]

scaf128 <- mapping_table[mapping_table$ref == "lg8_ord15_scaf128",]

scaf128[scaf128$block_r_start < 5.69e6 & scaf128$block_r_end > 5.16e6, ]
### generate blank variant_density_table
```

### Elevated heterozygosity in _T. californicum_

In following script I have done some exploration

```
E_colour_loci/estimate_Tcm_lg_8_color_loci_size.R
```

The length of the elevated heteoryzogisty is 23.8 Mbp.
