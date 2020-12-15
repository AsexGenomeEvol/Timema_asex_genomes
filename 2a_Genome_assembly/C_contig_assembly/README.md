I found that different assemblers perform better on different species. I chose three approaches that generated good assemblies for some species and tested them on all. The preselection of assembly approaches was done to reduce the computational load.

The three batches are :
- `SOAPdenovo2` - without additional sequencing, single-end and pair-end reads trimmed from mate-pair library.
- `SOAPdenovo2` - without single-end and pair-end reads trimmed from mate-pair library.
- `ABySS` - without additional sequencing, single-end and pair-end reads trimmed from mate-pair library.

The kmer choice for `SOAPdenovo2` batch1 was experimentally determined to be between 41 and 49. Therefore the optimal assembly was done by iterative search for local maxima in this range.

To assemble all the species with `SOAPdenovo2` and kmers 41 43 and 45 run (in the `2a_Genome_assembly` directory):

```
make assemble.batch1
```

`make asm.stats` will scan for all assemblies in `data/*/assemblies` and calculate formatted statics, `make asm.tables` will scan for stats that are not in the table of assemblies (one for every species in `stats/assemblies/<sp>_ctgs.tsv`). `make assemble.batch1.next_iter` will then look if assembly have reached local maxima and if not, calculate next kmer value.

```
make asm.stats # wait till finished, should not take more than couple of minutes
make asm.tables # very fast
make assemble.batch1.next_iter
```

The same principle with batch2 with only one difference, instead of initial exploration the first assembly take optimal values of batch1. Do not forget to recalculate stat before every next iteration, if tables wont get updated, iterative scripts will try to submit very same as in the last iteration (in fact it should report that assembly exists already)

```
make asm.stats
make asm.tables
make assemble.batch2.next_iter
```

## Testing of individual assemblers

I tested four species to access the optimal assembler for building contigs. Latter on I found that different assemblers perform better on different species (even they are just several million years divergent), therefore these test could not really find an optimal solution, they could only reduce parametric space for which I need to compute assembly of all the species. I ended up with three assembly approaches I call batch1, batch2 and batch3.

#### Naming conventions

The very basic stats I compute using my script, to avoid any misunderstanding from the reporting stats of assemblers and also to get unified format (which is totally retarded, my fault, but I wont change it for this project).

I have decided to make systematic naming of folders with assemblies, so I can see what the setting:

  category                          |   pattern in file name  |  meaning
 ---------------------------------- | ----------------------- | ------------------------------------------------------------------------
  software                          | `abyss`                 | assembly was performed by abyss
  .                                 | `SOAP`                  | assembly was performed by SOAPdenovo2
  .                                 | `megahit`               | assembly was performed by Megahit
  .                                 | `platanus`              | assembly was performed by Platanus
  k                                 | `_kXX`                  | kmer of assembly is `XX`
  kmer filter                       | (default)               | no threshold for filter was specified
  .                                 | `_kcX`                  | threshold for filter was specified to `X` (for abyss only)
  mate-pairs                        | (default)               | mate-pairs retrieved using NxTrim
  .                                 | `_fasteris`             | mate-pairs trimmed by fasteris    
  additional sequencing             | (default)               | input data include sequencing reads from December 2016
  .                                 | `_fewdata`              | input data include only original sequencing reads from February 2016
  single end                        | (default)               | all single end reads were used
  .                                 | `_nose`                 | no single end reads were used
  .                                 | `_nomse`                | only single end reads from pair end libraries were used (not those from mate pair library)
  pair reads from mate-pair library | (default)               | pair-end reads retrieved from mate-pair library were used
  .                                 | `_nompe`                | pair-end reads retrieved from mate-pair library were not used
  corrected reads                   |  (default)              | reads were not corrected
  .                                 | `_corrected`            | reads were corrected using `bfc`
  scaffolding method                | (default)               | native (for `scf` files) or none (for `ctg` files)
  .                                 | `_BESST`                | scaffolded using BESST
  .                                 | `_SOAP`                 | scaffolded using SOAPdenevo


Inside of every folder should be at least: folder with appropriate name and `*_ctg.stats` and/or `*_scf.stats` file - where these files are the output of the script `scripts/generic_genomics/fasta2genomic_stats.py`.

The other assemblers that got no systematic names were just computed on local disk and discarded afterwards.

#### SOAP

SOAPdenovo was used for the first assembly. Optimal kmer for SOAP was 41-49. The assembly requires a setting file, that can be created using script [[n3_setPaths.sh]] as follows: `C_contig_assembly/make_SOAP_setting.sh <sp> -nose -nompe -fewdata`.

```
SOAPdenovo-63mer all -s 2_Tcm.config -p 32 -a 150 -K 47 -R -o 2_Tcm_SOAP_k47 2> 2_Tcm_SOAP_k47.err
```

then the optimal kmergenie kmer was assembly was 83, so:

```
SOAPdenovo-127mer all -s 2_Tcm.config -p 32 -a 150 -K 83 -R -o 2_Tcm_SOAP_k83 2> 2_Tcm_SOAP_k83.err
```

The last thing to investigate is combined kmers:

```
SOAPdenovo-127mer all -s 2_Tcm.config -p 32 -a 150 -K 47 -m 83 -R -o 2_Tcm_SOAP_k47-83 2> 2_Tcm_SOAP_k47-83.err
```

This one was for some reason very very bad.

I wanted to see if SOAPdenovo produces a meaningful contigs even though it is using so much different kmers and it has so fragmented contigs compared to scaffolds. Therefore I gapfilled one assembly of T. californicum and evaluated it.

```
# dee-serv04
# SOAP asm with k=47
```

raw filtering before Gap Filling (to reduce computational time). For others we might want to do that mode elegantly (with proper remapping reads back).

```
grep ">" 2_Tcm_SOAP_k47_noncor.scafSeq | awk '{if($2 > 1){print $1}}' | sed s/.// > scf_to_keep.list
python3 $TROOT/scripts/generic_genomics/fasta2extract_by_list_of_headers.py 2_Tcm_SOAP_k47_noncor.scafSeq scf_to_keep.list > 2_Tcm_SOAP_k47_cov_fil.fa
python3 $TROOT/scripts/generic_genomics/fasta2genomic_stats.py 2_Tcm_SOAP_k47_cov_fil.fa > 2_Tcm_SOAP_k47_cov_fil.stats
python3 $TROOT/scripts/generic_genomics/fasta2fasta_length_filtering.py 2_Tcm_SOAP_k47_cov_fil.fa 500 > 2_Tcm_SOAP_k47_fil.fa
python3 $TROOT/scripts/generic_genomics/fasta2genomic_stats.py 2_Tcm_SOAP_k47_fil.fa 1300000000 > 2_Tcm_SOAP_k47_fil.stats
# Totally 231373 if 1397069692 nt. that's ok.
rm 2_Tcm_SOAP_k47_cov_fil.fa
```

and `Tcm_SOAP_k47_fil.fa` can be gap filled now (first locally to verify computational resources needed for the job).

```
/scratch/beegfs/monthly/ptranvan/Software/GapCloser/1.12-r6/GapCloser -a 2_Tcm_SOAP_k47_fil.fa -b 2_Tcm.config -o 2_Tcm_SOAP_k47_fil_GC.fa -t 32 -l 150
```

at one point it uses nearly 60G of memory. So, lets see how it worked:

```
module add UHTS/Quality_control/quast/4.1
quast.py 2_Tcm_SOAP_k47_fil_GC.fa -o 2_Tcm_SOAP_k47_fil_GC_quast -t 32 -s -f --eukaryote
```

took 8 hours. N looks ok (0.5% or something like that). I will run BUSCO as well.

```
BUSCO.py -i 2_Tcm_SOAP_k47_fil_GC.fa -o 2_Tcm_SOAP_k47_fil_GC_BUSCO -m geno -l /scratch/beegfs/monthly/kjaron/busco_ref/insecta_odb9 -c 32
```

Well, not that bad: C:85.0%[S:83.9%,D:1.1%],F:9.7%,M:5.3%,n:1658. There is quite a lot of genes fragmented (almost 10%), but the rest looks reasonable (not too much duplications, a lot of singletons).

###### More tests on Tdi

Tests in this subsectopm were done with repo in commit 9a6a9d

```bash
LOCALDIR=/scratch/local/monthly/kjaron/1_Tdi_SOAP_k43
mkdir -p $LOCALDIR/temp
export TMPDIR=$LOCALDIR/temp
cd $LOCALDIR
sed s/template/Tdi/g $TROOT/C_contig_assembly/template.config > Tdi.config
sed -i s:READPATH:$LOCALDIR:g Tdi.config
ln -s $TROOT/data/1_Tdi/trimmed_reads/*.fq.gz .
ln -s $TROOT/data/1_Tdi/trimmed_reads/mp_nxtrim_RF/*.fq.gz .
SOAPdenovo-63mer all -s Tdi.config -p 32 -a 120 -K 43 -R -o Tdi_SOAP_k43 &> 1_Tdi_SOAP_k43.err
```

continuity stats are worst then they were without more data. I will try to rerun the original assembly. And I found one more problem, that max read length was set to 125, instead of 150.

```
mkdir ../1_Tdi_SOAP_k43_fewdata
cd ../1_Tdi_SOAP_k43_fewdata
cp ../1_Tdi_SOAP_k43/Tdi.config ./Tdi_fewdata.config
```

manually delete 350_run2 libraries and trimmed se reads from mp library.

```
SOAPdenovo-63mer all -s Tdi_fewdata.config -p 32 -a 120 -K 43 -R -o Tdi_SOAP_k43_fewdata &> 1_Tdi_SOAP_k43_fewdata.err
```

In `Tcm` tests, I have not used any SE reads. SO I will try that as well

```
SOAPdenovo-63mer all -s Tdi_no_se.config -p 32 -a 120 -K 43 -R -o Tdi_SOAP_k43_no_se &> 1_Tdi_SOAP_k43_no_se.err
```

Both those assemblies (scaffolds) are less continuous than the raw assembly. Let's look at contig level. From the error log from the original raw assembly I diged:

```
There are 1807414 contig(s) longer than 100, sum up 1029303569 bp, with average length 569.
The longest length is 45913 bp, contig N50 is 1410 bp,contig N90 is 174 bp.
14039208 contig(s) longer than 44 output.
```

this is no single end, but with `350_run2` library:

```
There are 4674566 contig(s) longer than 100, sum up 1461866426 bp, with average length 312.
The longest length is 19883 bp, contig N50 is 614 bp,contig N90 is 142 bp.
18389576 contig(s) longer than 44 output.
```

and here is fewdata but with single end (only with pe):

```
There are 1688088 contig(s) longer than 100, sum up 1040643278 bp, with average length 616.
The longest length is 45913 bp, contig N50 is 1578 bp,contig N90 is 187 bp.
13046909 contig(s) longer than 44 output.
```

and this part is bit strange! Where the difference in contigs is comming from? The input data are very same (at least on contig level, mate pairs are different). However, it is a newer version of software, therefore small improvement might be expected (there is indeed bit less bit longer contigs).

Ok, now I can try very same combination of parameters for ABySS for kmer that is optimal for SOAP and kmer optimal for ABySS.


```bash
$TROOT/C_contig_assembly/2_SOAP.sh 1_Tdi 45 --nomse
$TROOT/C_contig_assembly/2_SOAP.sh 1_Tdi 45 --nomse --nompe
$TROOT/C_contig_assembly/2_SOAP.sh 1_Tdi 45 --nomse --nompe --fewdata
$TROOT/C_contig_assembly/2_SOAP.sh 1_Tdi 45 --nomse --nompe --fewdata --fasteris
$TROOT/C_contig_assembly/2_SOAP.sh 1_Tdi 45 --nomse --fasteris
```

#### ALLPATHS-LG

commands in the script [n1_allpaths.sh](n1_allpaths.sh)

program requires an overlapping library, which I have created from pair-end reads in mate-pair library (library `is_225`, see [B_read_parsing](../B_read_parsing)).
This library has too few overlaps to be used, therefore ALLPATHS-LG is disqualified from other testing.
Interesting: it evaluates a genome size, proportion of repetitive regions, etc...
Maybe it would be useful to use these modules for getting stats (those estimates).

#### Platanus

Command in the script [n2_Tps_platanus.sh](n2_Tps_platanus.sh). The assembly had in order of magnitude more assembled contigs than any other approach (~ 10M contigs, N50 238 bp, total sum 1.764 Gbp).

#### ABySS

Test on `Tps` (manually linked reads to local disk where assembly was performed). I also catted all R1 and R2 single end reads (which I have not done for the final assembly).

```sh
abyss-pe name=Tps_abyss_k83 k=83 lib='pea peb pec ped pee pef' pea='1_Tps_R1t_is_225.fq Tps_R2t_is_225.fq' peb='1_Tps_R1t_is_350.fq Tps_R2t_is_350.fq' pec='1_Tps_R1t_is_550.fq Tps_R2t_is_550.fq' ped='1_Tps_R1t_is_700.fq Tps_R2t_is_700.fq' pee='1_Tps_R1t_is_3000.fq Tps_R2t_is_3000.fq' pef='1_Tps_R1t_is_5000.fq Tps_R2t_is_5000.fq' se='1_Tps_se.fq'
```

since tests are showing in order of magnitude greater continuity of ABySS assembly than of `Platanus`,

Further testing of abyss was dome using bloom filter (it is [way faster and results are comparable](http://biorxiv.org/content/early/2016/08/07/068338))

```bash
mkdir -p /scratch/local/kjaron/temp
export TMP_DIR=/scratch/local/kjaron/temp
abyss-pe name=5_Tge_abyss_k83 j=32 k=83 lib='pea peb pec ped pef' mp='mpa mpb' \
 pea='5_Tge_R1t_is_225.fq.gz 5_Tge_R2t_is_225.fq.gz' \
 peb='5_Tge_R1t_is_350.fq.gz 5_Tge_R2t_is_350.fq.gz' \
 pec='5_Tge_R1t_is_350_run2.fq.gz 5_Tge_R2t_is_350_run2.fq.gz' \
 ped='5_Tge_R1t_is_550.fq.gz 5_Tge_R2t_is_550.fq.gz' \
 pef='5_Tge_R1t_is_700.fq.gz 5_Tge_R2t_is_700.fq.gz' \
 mpa='5_Tge_R1t_is_3000.fq.gz 5_Tge_R2t_is_3000.fq.gz' \
 mpb='5_Tge_R1t_is_5000.fq.gz 5_Tge_R2t_is_5000.fq.gz' \
 se='5_Tge_R1np.fq.gz 5_Tge_R2np.fq.gz' \
 B=26G H=4 kc=3 v=-v 1> abyss_corrected_Tcm_k83.log 2> abyss_corrected_Tcm_k83.err
```

will create assembly of `5_Tge` with kmer 83.

test with corrected reads:

```bash
mkdir -p /scratch/local/kjaron/temp
export TMP_DIR=/scratch/local/kjaron/temp
abyss-pe name=2_Tcm_abyss_k83 j=32 k=83 lib='pea peb pec ped pef' mp='mpa mpb' \
 pea='2_Tcm_R1tc_is_225.fq.gz 2_Tcm_R2tc_is_225.fq.gz' \
 peb='2_Tcm_R1tc_is_350.fq.gz 2_Tcm_R2tc_is_350.fq.gz' \
 pec='2_Tcm_R1tc_is_350_run2.fq.gz 2_Tcm_R2tc_is_350_run2.fq.gz' \
 ped='2_Tcm_R1tc_is_550.fq.gz 2_Tcm_R2tc_is_550.fq.gz' \
 pef='2_Tcm_R1tc_is_700.fq.gz 2_Tcm_R2tc_is_700.fq.gz' \
 mpa='2_Tcm_R1tc_is_3000.fq.gz 2_Tcm_R2tc_is_3000.fq.gz' \
 mpb='2_Tcm_R1tc_is_5000.fq.gz 2_Tcm_R2tc_is_5000.fq.gz' \
 se='2_Tcm_R1npc.fq.gz 2_Tcm_R2npc.fq.gz 2_Tcm_se_c_mp.fq.gz' \
 B=26G H=4 kc=3 v=-v 1> abyss_corrected_Tcm_k83.log 2> abyss_corrected_Tcm_k83.err
```

Note that these script I have written for testing were deleted. However, if you look at [this](https://github.com/AsexGenomeEvol/timema_assembly/tree/4fbf966671baf1059b7fd854f4a8ac9fb8369a36/C_contig_assembly) point of history, everything is there.

according conversation with developers of ABySS I should give a try multiple kmers and kc values. The script `1_abyss_grid.sh` should do the job. It is true that according kmerprofiles it is not clear at all what is the optimal kmer. In fact I guess that BWISE could do a good job, since it takes more kmers into account to construct superreads. To get an overview

```
cd $TROOT/data/2_Tcm/assembly
grep -A 2 -B 1 "unitigs.fa" *out
```

to get bit better overview:

```
cd $TROOT/data/2_Tcm/assembly
find . -name "*ctg*fa" -exec get_stats.sh {} \;
```

The stats collected by this, were parsed in script `../stats/ABySS_evaluation.R`. I have just spotted that I run the abyss in all the cases with non-paired reads, but all SOAP assemblies without. What more, the 225 library was also absent in SOAP assemblies. Damn! Assemblies without se reads and pair end mate pairs I will call "nose" (no single end).

```
$TROOT/C_contig_assembly/1_abyss_no_se.sh 2_Tcm 83 3
```

and test for `Tdi` comparable with SOAP denovo:

```
abyss-pe name=Tdi_abyss_k83_kc3_nose j=32 k=83 lib='peb ped pef' mp='mpa mpb'  peb='Tdi_R1t_is_350.fq.gz Tdi_R2t_is_350.fq.gz' ped='Tdi_R1t_is_550.fq.gz Tdi_R2t_is_550.fq.gz' pef='Tdi_R1t_is_700.fq.gz Tdi_R2t_is_700.fq.gz'  mpa='Tdi_R1t_is_3000.fq.gz Tdi_R2t_is_3000.fq.gz'  mpb='Tdi_R1t_is_5000.fq.gz Tdi_R2t_is_5000.fq.gz'  B=26G H=4 kc=3 v=-v 1> 1_Tdi_abyss_k87_kc3_fewdata_nose_nompe.out 2> 1_Tdi_abyss_k87_kc3_fewdata_nose_nompe.err
```

#### Megahit

should be a successor of SOAP, therefore I will try it even though it is designed for metagenomic assemblies. Min count should have similar properties as kc in ABySS, therefore I know that 3 seems the optimal value at least for 2_Tcm.

```bash
cd /scratch/local/kjaron/2_Tcm_corrected_reads/Tcm_megahit_k47-83_run1
megahit -1 2_Tcm_R1t_is_225.fq.gz -2 2_Tcm_R2t_is_225.fq.gz \
        -1 2_Tcm_R1t_is_350_run2.fq.gz -2 2_Tcm_R2t_is_350_run2.fq.gz \
        -1 2_Tcm_R1t_is_350.fq.gz -2 2_Tcm_R2t_is_350.fq.gz \
        -1 2_Tcm_R1t_is_550.fq.gz -2 2_Tcm_R2t_is_550.fq.gz \
        -1 2_Tcm_R1t_is_700.fq.gz -2 2_Tcm_R2t_is_700.fq.gz \
        -r 2_Tcm_R1np.fq.gz -r 2_Tcm_R2np.fq.gz -r 2_Tcm_se_mp.fq.gz \
        --min-count 3 --k-min 47 --k-max 83 --k-step 6 \
        -m 0.2 -t 32 --out-prefix 2_Tcm_megahit_k47-83 \
        --tmp-dir /scratch/local/kjaron/temp --verbose \
        2> 2_Tcm_megahit_k47-83.err
```

so far it is using 60G of memory. It generated  file `megahit_out/Tcm_megahit_k47-83.contigs.fa`. Stats are not that bad as SOAP stats, but not that good that ABySS stats. The total length is bit smaller than expected (1.2G), but still reasonable. Maybe it could be fixed by using smaller filter for noisy kmers (2 instead of 3). Now finally I can just scaffold it (in section D).

The second try was done without 225 and single end libraries. Also step was chosen 12 instead of original 6.

#### Some notes

When contigs are not filtered, the scaffolding is super crazy (discussed in section D). I will filter on individual basis contigs shorter than certain length. The threshold will probably be 200 in the most of the cases. This step is usually done internally in assemblers that contain module for scaffolding.
