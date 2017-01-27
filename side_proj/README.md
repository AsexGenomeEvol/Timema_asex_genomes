### Variant calling: SNPs & Indels

so far I was thinking of FreeBayes for variant calling, but according the latest benchmarking I should probably go for Octopus.

Mapping using `bwa-mem`

```
module add UHTS/Aligner/bwa/0.7.13

cd $TROOT/data/5_Tge/reference

# only scaffolds longer than 100k are kept (to make it computable)
python3 $TROOT/scripts/generic_genomics/fasta2fasta_length_filtering.py Tge_abyss87_besst_GC_core_scaffolds.fa 100000 > Tge_abyss87_besst_GC_core_100k_filtered.fa

# the total size of this reference is ~600M.
bwa index Tge_abyss87_besst_GC_core_100k_filtered.fa

$TROOT/scripts/map_pair_end_RG_lsf.sh 5_Tge \
  Tge_abyss87_besst_GC_core_100k_filtered.fa \
  Tge_R1t_is_350.fq.gz \
  Tge_R2t_is_350.fq.gz \
  "@RG\tID:TGE\tSM:TGE_REF\tPL:illumina\tLB:350\tPU:lane1"
```

and similar for `1_Tps`:

```
cd $TROOT/data/1_Tps/raw_assembly

# only scaffolds longer than 100k are kept (to make it computable)
python3 $TROOT/scripts/generic_genomics/fasta2fasta_length_filtering.py 1_Tps_genome.fa 100000 > 1_Tps_genome_100k_filtered.fa

# the total size of this reference is ~200M, 1326 scaffolds.
bwa index 1_Tps_genome_100k_filtered.fa

$TROOT/scripts/map_pair_end_RG_lsf.sh 1_Tps \
  Tge_abyss87_besst_GC_core_100k_filtered.fa \
  Tps_R1t_is_350.fq.gz \
  Tps_R2t_is_350.fq.gz \
  "@RG\tID:TPS\tSM:TPS_REF\tPL:illumina\tLB:350\tPU:lane2"
```

I am supposed to mark duplicates, but I wont, seems like a pint I do not have time for it now. However, script `1_deduplicate_lsf.sh` should do it (not tested).

Now we have to index `.fasta` and run [GATK](https://software.broadinstitute.org/gatk) (as a quick solution for now).
```
cd $TROOT/data/1_Tps/reference
$TROOT/scripts/index_fa.sh 1_Tps_genome_100k_filtered.fa
$TROOT/scripts/make_dict_fasta.sh 1_Tps_genome_100k_filtered.fa
cd ../variant_calling/GATK/350/
$TROOT/scripts/index_bam.sh map_pe_to_1_Tps.bam
# once those jobs are done
$TROOT/N_variant_calling/2_run_GATK_lsf.sh 1_Tps 350 1_Tps_genome_100k_filtered.fa
```
analogically I run the analysis on Tge, both 350 and 550 libs. For tps it finished without producing any result, it is not very clear why. I do not really like java programs so I try [atlas](https://bitbucket.org/phaentu/atlas/wiki/Home).

Seems, that estimates of theta are pretty light weighted. However, I should really sit and read how do they do that. Simply:

```
atlas task=callBayes bam=../../GATK/350/map_pe_to_1_Tps.bam 1> callBayes.out 2> callBayes.err
```
