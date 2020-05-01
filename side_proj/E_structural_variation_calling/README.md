## Structural variations

Detection of structural variants from population data will be based on pair-end illumina data.
One individual per species will have mate-pairs as well to access structural variations
with greater resolution.

### Inference from Illumina reads

Methods are based on coverage or mapping patterns. Different methods produce different variant calls. This observation let to development of methods that calls consensus based on different methods. The consensus caller is called (Survivor)[https://github.com/fritzsedlazeck/SURVIVOR].

#### SV calls: pair-end


Restarting efforts at Edinburgh required quite lot of logistic investment to move the data here. I have them in agregated directories by type rather than indevidual as before. Not sure which is better. If I will reorganize them again, I will have rewrite this, but for now I keep it as it is.

**data**

- 60 mapping files: `data/mapped_reseq_reads/*.bam`
- 60 manta SV calls (called at Vital-it cluster using `E_structural_variation_calling/manta.sh`): `data/manta_SV_calls/<sp>_<ind>_manta/*.vcf.gz`
- 60 lumpy SV calls (called at Vital-it cluster using `E_structural_variation_calling/smoove.sh`; smoove is a wrapper round lumpy): `data/smoove_SVs/...`
- reference genomes: `data/final_references/*_b3v08.fasta.gz`

```
for bam in data/mapped_reseq_reads/*.bam; do
  qsub -o logs/ -e logs/ -cwd -N bamindex -V -pe smp64 1 -b yes "samtools index $bam"
done
```

### Done in Lausanne

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

For a record, Manta does not by default call inversions, but breakpoints: `$MANTA_INSTALL_FOLDER/libexec/convertInversion.py` should convert BND to INV.

##### Delly

Check their [readme](https://github.com/dellytools/delly#germline-sv-calling), it's very clear and simple. installed it is running

```
delly call -g Tge_abyss87_besst_GC_core_100k_filtered.fa \
  map_pe_to_5_Tge.bam -o tge_350_delly.bcf
```

produced in 2 hours a `.bcf` file that was not parsed so far.

##### Lumpy

I installed lumpy and [smoove](https://brentp.github.io/post/smoove/) which is a probabilistic genotyper using SV calls by lumpy.

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

#### SV calls: mate-pairs

I have mate-pairs for only reference genomes.

#### Merging calls

Right now I have Delly, Lumpy and Manta SV calls, the union in T. monikensis is ~10k, overlap of at least two >4k and all three >1k. There are two strategies I will consider:

1. Accept all calls made by at least two callers
2. Create a union of all calls, and use a genotyper to test this set of candidate SV in all the individuals.

In either the case I will remove all calls homozygous in all (nearly all?) individuals (asm errors).

##### Making union

There two ways how to make a union. Delly or SURVIVOR. Delly uses both recoprocal overlap and breakpoint offset to consider an SV the same. SURVIVOR focuses on the offset only. Delly is however screwing up on merging other SV callers, so SURVIVOR it is

```
E_structural_variation_calling/survivor_all_merged_calls.sh <sp>
```

generates

```
data/$SP/variant_calls/"$SP"_survivor_all_calls_union.vcf
```

file with the default merging parameters (min len > 30; breakpoint distance < 1000). It is a wild script for now, but once I will have delly SV calls for all the species, I will embed it to `Snakemake`.

- individual delly / smove / manta SV calls (mostly if we want to go back and check something)
- merged calls of the three (`"$SP"_all_calls_merged.bcf`, is there support inside? Need to check)

#### Genotyping

This step is now essential. Till now we were calling variants with a crtain reliability. We need to take their union and ask again, in what samples is this variant present? Absence of clear evidence (the reason why it was not called in the first place) and we need to confirm the absence rather.

##### by Delly

Delly has a genotyper given set of candidates, so I use it while feeding it with the merged SVs.

```
E_structural_variation_calling/delly_genotyping.sh <sp>
```

generating

```
data/$SP/variant_calls/$SAMPLE/delly_genotyping.bcf # for each sample
data/$SP/variant_calls/delly_genotyping_merged.bcf
```

Output are merged genotyping calls (`"$SP"_delly_genotyping_merged.bcf`), given the set of candidates

### Done in Edinburgh

##### Paragraph

This is program of choice. The problem is that it requires formating for the `.vcf` file that contains `REF` and `ALT` (check minimalist example `data/testing_data/round-trip-genotyping/candidates.vcf`). Which means that I need to go one step back, to figure out how to merge calls WITH explicit `REF`/`ALT` sequences. The other option would be to genotype using calls of each other and then base the mergin on the genopyping itself (SURVIVOR on steroids). I should try to genotype reciprocally two samples and then see if the genotyping is consistent with SURVIVOR merged calls.

I feed `paragraph` with processed manta SV calls using `scripts/convertManta/convertManta2Paragraph_compatible_vcf.py` script.

```
parallel -j 1 'qsub -o logs/ -e logs/ -cwd -N paragraph -V -pe smp64 32 -b yes {}' :::: E_structural_variation_calling/paragraph_commands.txt
```

This works now for `INS`, `DEL`. Does not for `DUP` and `INV`. The full log of problems I have encountered is bellow.


##### Paragraph error log

- Error 1: [SVs too close to the start of scafflds](https://github.com/Illumina/paragraph/issues/42). For now I will just kick them out
- Error 2: [Unresloved <INS>](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md#insertions-with-incomplete-insert-sequence-assembly) by definition don't have resolved sequence. I will filter those out as well for now. I guess the better strategy would be to use Ns.
- Error 3: `Illegal character in ALT allele: [3_Tms_b3v08_scaf000202:631916[N` I think this means it does not like the the inversions!
- Error 4: `Exception: 514459:514542 missing REF or ALT sequence.`; and the problematic SV is `3_Tms_b3v08_scaf000253	514459	MantaINV:5899:0:0:0:0:0	N	<INV>	34	SampleFT	END=514542;SVTYPE=INV;SVLEN=83;IMPRECISE;CIPOS=-351,352;CIEND=-355,356;INV5	GT:FT:GQ:PL:PR	1/1:MinGQ:9:86,12,0:0,4`. Reported [here](https://github.com/Illumina/paragraph/issues/43). For now I just manually deleted it, as I don't see anything wrong with it for now.
- Error 5: `Exception: Illegal character in ALT allele: ]3_Tms_b3v08_scaf007303:7862]T` This is a breakpoint that is not an inversion, so it was not convered.

All the are now in a single script `scripts/convertManta/convertManta2Paragraph_compatible_vcf.py`:

```
conda activate default_genomics
SCRIPT=scripts/convertManta/convertManta2Paragraph_compatible_vcf.py
REF=data/final_references/3_Tms_b3v08.fasta.gz
VARIANTS=data/manta_SV_calls/Tms_00_manta/results/variants/diploidSV.vcf.gz
OUT=data/manta_SV_calls/Tms_00_manta/results/variants/diploidSV_corrected_decomposed
python3 $SCRIPT $REF $VARIANTS --quality_filtering -prefix_split_by_type $OUT
```

**INV**

The current warnings:

```
2020-04-06 15:54:52,957 WARNING  3_Tms_b3v08_scaf000046:667050 Padding base in genome is different from VCF. Use the one from genome.
2020-04-06 15:54:53,125 WARNING  3_Tms_b3v08_scaf000820:160049 Padding base in genome is different from VCF. Use the one from genome.
2020-04-06 15:54:53,128 WARNING  3_Tms_b3v08_scaf000273:75266 Padding base in genome is different from VCF. Use the one from genome.
...
2020-04-06 15:54:58,956 WARNING  3_Tms_b3v08_scaf009637:199 Padding base in genome is different from VCF. Use the one from genome.
```

I should figure out why the warnings and potentially fix the conversion script.

31/72 variants get the warning.

```
3_Tms_b3v08_scaf000820  160049  MantaINV:12347:0:0:1:0:0        B'>3_TMS_B3V08_SCAF000820:160049-160049\NN\N'   <INV>   28      PASS    END=160307;SVTYPE=INV;SVLEN=258;IMPRECISE;
3_Tms_b3v08_scaf000273  75266   MantaINV:6190:0:0:2:0:0 B'>3_TMS_B3V08_SCAF000273:75266-75266\NC\N'     <INV>   87      PASS    END=76145;SVTYPE=INV;SVLEN=879;IMPRECISE;CIPOS=-312
...
```

These are only those that are marked as `IMPRECISE`.

Compared to those that are just alright:

```
3_Tms_b3v08_scaf000266  245751  MantaINV:6090:0:0:0:0:0 T       <INV>   32      PASS    END=245911;SVTYPE=INV;SVLEN=160;IMPRECISE;CIPOS=-326,327;CIEND=-375,375;INV3
```

Examples of two types of converted inversions:

##### GOOD

```
3_Tms_b3v08_scaf000034  455072  MantaBND:1196:0:1:0:0:0:1       T       T]3_Tms_b3v08_scaf000034:455954]        29      PASS    SVTYPE=BND;MATEID=MantaBND:1196:0:1:0:0:0:0;IMPRECI
3_Tms_b3v08_scaf000034  455954  MantaBND:1196:0:1:0:0:0:0       G       G]3_Tms_b3v08_scaf000034:455072]        29      PASS    SVTYPE=BND;MATEID=MantaBND:1196:0:1:0:0:0:1;IMPRECI
```

into

```
3_Tms_b3v08_scaf000034  455072  MantaINV:1196:0:1:0:0:0 T       <INV>   29      PASS    END=455954;SVTYPE=INV;SVLEN=882;IMPRECISE;CIPOS=-318,318;CIEND=-319,319;INV3    GT:FT:GQ:PL
```

4. 3_Tms_b3v08_scaf000076

```
X  206212  MantaBND:2416:0:1:0:0:0:0       A       A]X:209142]        122     PASS    SVTYPE=BND;MATEID=MantaBND:2416:0:1:0:0:0:1;CIPOS=0
X  209131  MantaBND:2416:0:1:0:0:0:1       T       T]X:206223]        122     PASS    SVTYPE=BND;MATEID=MantaBND:2416:0:1:0:0:0:0;CIPOS=0
```

into

```
X  206212  MantaINV:2416:0:1:0:0:0 A       <INV>   122     PASS    END=209142;SVTYPE=INV;SVLEN=2930;CIPOS=0,11;CIEND=-11,0;HOMLEN=11;HOMSEQ=ATAGAGGTTAG;INV3
```

#### BAD

```
3_Tms_b3v08_scaf000001  1881491 MantaBND:49:0:0:1:0:0:0 C       [3_Tms_b3v08_scaf000001:1881608[C       43      PASS    SVTYPE=BND;MATEID=MantaBND:49:0:0:1:0:0:1;CIPOS=0,2;HOMLEN=
3_Tms_b3v08_scaf000001  1881606 MantaBND:49:0:0:1:0:0:1 A       [3_Tms_b3v08_scaf000001:1881493[A       43      PASS    SVTYPE=BND;MATEID=MantaBND:49:0:0:1:0:0:0;CIPOS=0,2;HOMLEN=
```

into

```
3_Tms_b3v08_scaf000001  1881490 MantaINV:49:0:0:1:0:0   B'>3_TMS_B3V08_SCAF000001:1881490-1881490\NG\N' <INV>   43      PASS    END=1881607;SVTYPE=INV;SVLEN=117;CIPOS=0,2;CIEND=-2
```

**3_Tms_b3v08_scaf000046**

```
X  667051  MantaBND:1568:0:0:1:0:0:1       G       [X:667176[G        40      PASS    SVTYPE=BND;MATEID=MantaBND:1568:0:0:1:0:0:0;IMPRECI
X  667176  MantaBND:1568:0:0:1:0:0:0       A       [X:667051[A        40      PASS    SVTYPE=BND;MATEID=MantaBND:1568:0:0:1:0:0:1;IMPRECI
```

into

```
X  667050  MantaINV:1568:0:0:1:0:0 B'>3_TMS_B3V08_SCAF000046:667050-667050\NT\N'   <INV>   40      PASS    END=667175;SVTYPE=INV;SVLEN=125;IMPRECISE;CIPOS=-38
```

This seems like the problem happens then there is reverse complementary reference sequence. Let's see why is that.

Some other breakpoints entries

```
3_Tms_b3v08_scaf000018  369206  MantaBND:658:0:0:2:0:0:0        A       [3_Tms_b3v08_scaf000018:369560[AAAACAAACAAGTAAAAACAAGTA 698     PASS    SVTYPE=BND;MATEID=MantaBND:658:0:0:
3_Tms_b3v08_scaf000018  369560  MantaBND:658:0:0:2:0:0:1        G       [3_Tms_b3v08_scaf000018:369206[ACTTGTTTTTACTTGTTTGTTTTG 698     PASS    SVTYPE=BND;MATEID=MantaBND:658:0:0:
3_Tms_b3v08_scaf000019  938566  MantaBND:707:0:0:0:0:0:1        G       G]3_Tms_b3v08_scaf000019:938717]        15      MinQUAL SVTYPE=BND;MATEID=MantaBND:707:0:0:0:0:0:0;IMPRECIS
3_Tms_b3v08_scaf000019  938717  MantaBND:707:0:0:0:0:0:0        G       G]3_Tms_b3v08_scaf000019:938566]        15      MinQUAL
```

Questions:
 1. Why the `3_Tms_b3v08_scaf000018` inversion got filtered? Is it an inversion?

Ha, that's because it's one of the inversions with `SVINSSEQ=AAAACAAACAAGTAAAAACAAGT`. That means that on one side of the inversion there is an insertion too. I am for now filtering those variants as paragraph does not have a capacity to genotype nested variants.


**DUP**

Paragraph README:

```
<DUP> for duplication
    Must have END key in INFO field. paraGRAPH assumes the sequence between POS and END being duplicated for one more time in the alternative allele.
```

too many erros. I will try to do it locally now:

here I am testing only `PASS` duplications

```
OUT=data/manta_SV_calls/Tms_00_manta/results/variants/diploidSV_corrected_decomposed_quality_filtered
python3 $SCRIPT $REF $VARIANTS -prefix_split_by_type --quality_filtering $OUT
multigrmpy.py -i data/manta_SV_calls/Tms_00_manta/results/variants/diploidSV_corrected_decomposed_DUP.vcf -m data/genotyping/Tms_samples.txt -r data/final_references/3_Tms_b3v08.fasta.gz -o data/genotyping/3_Tms_ind00_genotyping_DUP --scratch-dir temp/3_Tms_ind00_genotyping_DUP --threads 16
```

**other**

Samples requires depth of each bam file, so

```
samtools depth -a data/mapped_reseq_reads/Tms_00_to_b3v08_mapped_within_scfs.bam > coverage_depth # coule be directly piped
cat coverage_depth |  awk '{sum+=$3} END { print "Average = ", sum/NR}'
## 11.55
## av read len?? made it up to 100
```

Get coverages of all bam files

```
qsub -o logs/ -e logs/ -cwd -N coverages -V -pe smp64 1 -b yes "bash calculate_coverages.sh"
```

##### GraphTyper2



#### TO CONSIDER

I think smove and manta use different names for the same thing (duplication vs insertion) or at least they sums are the same and one distinguishes them and the other does not. So it might be a good idea to "unify" them before merging. SURVIVOR cared about SV typpes, not sure how exactly Delly merger works.

- [stix](https://github.com/ryanlayer/stix)

CNV with a different specialised software, such as http://software.broadinstitute.org/software/genomestrip/
