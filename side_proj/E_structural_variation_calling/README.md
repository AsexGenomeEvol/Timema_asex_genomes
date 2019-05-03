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