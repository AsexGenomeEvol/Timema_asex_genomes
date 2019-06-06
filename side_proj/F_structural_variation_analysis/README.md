## Analysis of SVs

Analysis of 6 individuals per species. I want to know heterozygosity patterns in individual species. Here we do:

 - basic population genetics on SVs - site frequency spectra, heterozygosity analysis
 - the density of SVs on scaffolds (addressing putative reproduction mode)
 - per type / per size analysis of the SVs in the genomes.

#### What do we have in the end

- individual delly / smove / manta SV calls (mostly if we want to go back and check something)
- merged calls of the three (`"$SP"_all_calls_merged.bcf`, is there support inside? Need to check)
- merged genotyping calls (`"$SP"_delly_genotyping_merged.bcf`), given the set of candidates

### basic population genetics

- SFS (site freq. spectra)
- SFS vs heterozygosity (figure out how to plot it, maybe monochromatic scale and a simple tile per type??)
- stats as fractions of "present in all", "heterozygous in all", "homozygous in all", "heterozygous"

lot of this is done in very unsorted way in `playground.R` script.

### density of SVs

- check modality in sexual / asexuals (sexual are expected to segregate regardless of position in the genome, i.e. all scaffolds should be even)
- in asexuals we expect unimodal disr for homozygous SVs and bimodal for heterozygous SVs. Furthermore, we expect a negative relation of the two
- A statistical approach? Given a random placement, what is the expected distribution of observed density given the scf size, quantile?

### per Type / Size analysis

- For this I should probably separate homozygous and heterozygous SVs in asexuals (I will have 5 inds for homo, 6 for hetero).