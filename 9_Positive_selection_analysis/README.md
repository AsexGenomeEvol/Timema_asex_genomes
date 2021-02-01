# Positive_selection


We ran the [Selectome](https://selectome.org/timema) positive selection pipeline on the 7157 orthologs identifed in [3_Orthologs](../3_Orthologs). 

Full methods are available on the selectome website, but breifly:

All alignment building and filtering was performed on predicted amino acid sequences, and the final amino acid MSAs (multiple sequence alignments) were used to infer the nucleotide MSAs used for positive selection inference. MSAs were obtained by MAFFT (v. 7.310) with the allowshift option, which avoids over-aligning non homologous regions (e.g. gene prediction errors, or alternative transcripts). All the next steps "mask" rather than remove sites, by replacing the amino acid with a 'X' and the corresponding codon with ‘NNN’. MCoffee (v11.00.8cbe486)  was run with the following aligners: mafft_msa, muscle_msa, clustalo_msa, and t_coffee_msa. MCoffee provides a consistency score per amino acid, indicating how robust the alignment is at that position for that sequence. Residues with a consistency score less than 5 were masked. TrimAl (v. 1.4.1) was used to mask columns with less than 4 residues (neither gap nor 'X'). Following this step 2 of the 7157 ortholog alignments consisted only of gaps and were excluded from further analyses. 

Alignments available [here](https://selectome.org/timema/download).

The branch-site model with rate variation at the DNA level was run using the Godon software (https://bitbucket.org/Davydov/godon/, version 2020-02-17, option BSG --ncat 4). Each branch was tested iteratively, in one run per gene tree. For each branch, we obtain a ΔlnL which measures the evidence for positive selection, a corresponding p-value and associated q-value (estimated from the distribution of p-values over all branches of all genes), and an estimate of the proportion of sites under positive selection if any.



