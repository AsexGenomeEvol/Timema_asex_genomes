## Structural variations

Detection of structural variants from population data will be based on pair-end illumina data.
One individual per species will have mate-pairs as well to access structural variations
with greater resolution.

### Inference from Illumina reads

Methods are based on coverage or mapping patterns. Different methods produce different variant calls. This observation let to development of methods that calls consensus based on different methods. The consensus caller is called (Survivor)[https://github.com/fritzsedlazeck/SURVIVOR].

#### pair-end

##### Manta

The tools is designed for tumor cells, but there is no reason to think that it wont work on whole bodies.

I downloaded precompiled binaries of (Manta)[https://github.com/Illumina/manta] v1.0.3.

Unpacking and running `runMantaWorkflowDemo.py` reported succesful run and some details of commends that were tested. I will extract those commands and apply them to Tge, since I already have indexed reference and mapped pair end reads..

```
/Home/kjaron/src/manta-1.0.3.centos5_x86_64/bin/configManta.py --normalBam="map_pe_to_5_Tge.bam"  --referenceFasta="Tge_abyss87_besst_GC_core_100k_filtered.fa" --runDir=$(pwd)
```

the workflow script was created, so I run it. It can nicely limit the memory.

```
python -E runWorkflow.py -m local -j 32 -g 120 1> tge_350_manta.log 2> tge_350_manta.err
```

53k SVs, 7k of good quality. Very small Indels included.

Turned into scripts:

The name of the reference individual's insert size 350: `ref_is350`

```bash
manta.sh 1_Tdi b3v04 ref_is350
```

##### Delly

Check their [readme](https://github.com/dellytools/delly#germline-sv-calling), it's very clear and simple. installed it is running

```
delly call -g Tge_abyss87_besst_GC_core_100k_filtered.fa \
  map_pe_to_5_Tge.bam -o tge_350_delly.bcf
```

produced in 2 hours a `.bcf` file that was not parsed so far.

##### Lumpy

I installed lumpy and [smoove](https://brentp.github.io/post/smoove/) which is a probabilistic genotyper using SV calls by lumpy.

```
module add UHTS/Analysis/samtools/1.8
```

test run:

```
fasta=data/2_Tcm/reference/2_Tcm_b3v08.fasta.gz
bam=data/2_Tcm/mapping/Tcm_05_to_b3v08.bam
outdir=data/2_Tcm/variant_calls/Tcm_05/Tcm_05_smoove/
smoove call -x --genotype --name Tcm_05 --outdir $outdir \
           -f $fasta --processes 24 $bam
```


##### Brakedancer

weird problem with libraries.

#### mate-pairs

I have mate-pairs for only reference genomes

#### Merging calls

Right now I have Delly, Lumpy and Manta SV calls, the union in T. monikensis is ~10k, overlap of at least two >4k and all three >1k. There are two strategies I will consider:

1. Accept all calls made by at least two callers
2. Create a union of all calls, and use a genotyper to test this set of candidate SV in all the individuals.

In either the case I will remove all calls homozygous in all (nearly all?) individuals (asm errors).

#### Making union

There two ways how to make a union. Delly or SURVIVOR. Delly uses both recoprocal overlap and breakpoint offset to consider an SV the same. SURVIVOR focuses on the offset only. Might be a good idea to try both, as it brings only very little effort.

```
E_structural_variation_calling/delly_all_merged_calls.sh <sp>
```

generates

```
data/$SP/variant_calls/all_calls_merged.bcf
```

file with the default merging parameters (covergage > 10; overlap > 80%; max offset < 1000). It is a wild script for now, but once I will have delly SV calls for all the species, I will embed it to `Snakemake`.

#### Genotyping by Delly

Delly has a genotyper given set of candidates, so I use it while feeding it with the merged SVs.

```
E_structural_variation_calling/delly_genotyping.sh <sp>
```

generating

```
data/$SP/variant_calls/$SAMPLE/delly_genotyping.bcf # for each sample
data/$SP/variant_calls/delly_genotyping_merged.bcf
```

#### Paragraph

TODO

#### What do we have in the end

- individual delly / smove / manta SV calls (mostly if we want to go back and check something)
- merged calls of the three (`"$SP"_all_calls_merged.bcf`, is there support inside? Need to check)
- merged genotyping calls (`"$SP"_delly_genotyping_merged.bcf`), given the set of candidates

#### TO CONSIDER

I think smove and manta use different names for the same thing (duplication vs insertion) or at least they sums are the same and one distinguishes them and the other does not. So it might be a good idea to "unify" them before merging. SURVIVOR cared about SV typpes, not sure how exactly Delly merger works.