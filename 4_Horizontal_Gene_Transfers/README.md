# Horizontal gene transfers

Content:
  * **contamination_sequences** | Suspected contamination sequences identifed via blobtools (removed from the reference assemeblies prior to submission to NCBI).


To detect HGT from non-metazoan species, we first used the pipeline of foreign sequence detection developed by [Francois et al. 2020](https://www.g3journal.org/content/10/2/721). We used the set of CDS identified in publicly available transcriptomes and the genome assemblies prior to the decontamination procedure with Blobtools (see [2a_Genome_assembly](../2a_Genome_assembly) section). The rationale is that some genuine HGT could have been wrongly considered as contaminant sequences during this decontamination step and thus been removed from the assembly. Scaffolds filtered during decontamination are available in directory [contamination_sequences](contamination_sequences).

The pipeline is run via two bash scripts (which call the accessory scripts):

1. **script1_similarity_synteny_clustering_completion.sh**
* performs for each Timema species (independently) the similarity and synteny steps -> HGT candidates
* clustering of all HGT candidates identified in the 10 Timema species -> putative HGT families
* completion of all HGT families ('rescue' of homologs not yet identified in the genomic scaffolds) -> completed putative HGT families
* retrieve the reference ('blast-hitting') sequences for each HGT family -> ready for alignment.


2. **script2_phylogenetic_validation.sh** starts from the completed putative HGT families and performs the automated alignment then the tree reconstruction using RAxML.


Briefly, a DIAMOND BlastP (v0.8.33) allows to detect candidate non-metazoan genes in the set of CDS of each species. Taxonomic assignment is based on the 10 best blast hits to account for potential contaminations and other sources of taxonomic misassignment in the reference database. Candidate non-metazoan sequences are then subjected to a synteny-based screen with Gmap (v2016-11-07) to discriminate between contaminant sequences and potential HGT-derived sequences. A sequence is considered as a HGT candidate if it is physically linked to (i.e., mapped to the same scaffold as) at least one “confident-arthropod” CDS (previously identified in the DIAMOND  blast).
We then clustered all HGT candidates identified in each of the 10 Timema species into HGT families using Silix (v1.2.10), requiring a minimum of 85% identity (default parameters otherwise). These HGT families were then “completed” as much as possible by adding homologs from the genome assemblies not identified as HGT candidates (this could occur if the corresponding sequences are fragmented or on short scaffolds for example). To this end, the longest sequence of each HGT family was mapped (using Gmap) on the genomic scaffolds of all species, requiring a minimum of 85% identity.
For each completed HGT family, a protein alignment of the candidate HGT sequence(s) and its (their) 50 best DIAMOND blastP hits in the reference database (1st step of the pipeline) was generated with MAFFT (v7). The alignments were cleaned using HMMcleaner (stringency parameter = 12) and sites with more than 50% missing data were removed. Phylogenetic trees were inferred using RAxML (v8.2) with the model ‘PROTGAMMALGX’ of amino-acid substitution and 100 bootstrap replicates. Phylogenetic trees were inspected by eye to confirm or not an evolutionary history consistent with the hypothesis of HGT.



