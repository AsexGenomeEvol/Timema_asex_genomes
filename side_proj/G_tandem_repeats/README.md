### Tandem repetitions

Inspired by this [paper](https://academic.oup.com/mbe/article/35/4/925/4817473). Tandem repetitions can evolve rapidly. They also very apparently related to recombination as the mechanism that changes the number of copies between haplotypes (better explained in the paper).

### Annotation

There are probably tons of methods to annotate tandem repeats. However, we use the one they have used in the paper. It does worry me, as it does not seem to be picked up by the community and therefore it's hard to say whether the number we will get are of any use.

However, we can do a sanity check. Running the software of all the reseq data will allow us to estimate variation in the method. We know that `5_Tge` individuals are nearly clones and therefore we can expect that the differences in annotations will be only due to methodology. Hopefully we will be able to get about of changes in tandem repeats across the _Timema_ tree.

This trial have worked:

```
bsub -M 50000000 -q bgee -o kseek_trial.log 'perl /home/kjaron/bin/k_seek.pl Tdi_550.fastq 1_Tdi'
```

Now, I want to generate an annotation for all... This script

```
G_tandem_repeats/tandem_repeats.sh
```

and rules in snakefile do the job. Now I have `data/<sp>/tandem_repeats/<ind>_kseek.rep.total` files. The exploratory analysis will be done in `G_tandem_repeats/analysis_of_simple_repeats.R` script.