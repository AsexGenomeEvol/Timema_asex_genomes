# Orthologs

Description from the paper:

> Timema orthologous groups (OGs) were inferred with the OrthoDB standalone pipeline (v. 2.4.4) using default parameters. In short, genes are clustered with a graph-based approach based on all best reciprocal hits between each pair of genomes. The high level of fragmentation typical for Illumina-based genomes constrains the ability to identify 1:1 orthologs across all ten Timema species. To maximize the number of single copy OGs covering all ten Timema species, transcriptomes were included during orthology inference. Thus, transcripts were used to complete OGs in absence of a gene from the corresponding species. Using this approach, 7157 single copy OGs covering at least three sexual-parthenogenetic sister species pairs were obtained (SM Table 4).

# Input data

The input data was longest isoform amino acid sequences from both the genomes and previously published transcripomes.

### Genomes (input_data/genome_longest_isoforms.tar.gz)

Amino acid sequences were extracted from the genome assemblies using gffread (from Cufflinks v.2.2.1 (Trapnell et al. 2012)).

```
gffread sp.gff -g  sp_genome.fa  -y sp_aa.fa
```

### Transcriptomes (input_data/transcriptome_longest_isoforms.tar.gz)

Transcriptome sequences were taken from Bast et al. (2018). 



# TODO:
add a table of orthologs and link it


# References 

Bast, J., Parker, D. J., Dumas, Z., Jalvingh, K. M., Tran Van, P., Jaron, K. S., Figuet, E., Brandt, A., Galtier, N., Schwander, T. 2018. Consequences of asexuality in natural populations: insights from stick insects. Molecular Biology and Evolution. 35: 1668–1677.

Trapnell C, Roberts A, Goff L, Pertea G, Kim D, Kelley DR, et al. 2012. Differential gene and transcript expression analysis of RNA-seq experiments with TopHat and Cufflinks. Nat Protoc 7: 562–578.
