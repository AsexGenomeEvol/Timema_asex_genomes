#!/usr/bin/env bash

USE_LOCAL="yes" snakemake $@ -p --jobs 10 --cluster "bsub \
    -J {rule} \
    -q bgee \
    -n {threads} \
    -M {resources.mem} \
    -R \"span[hosts=1] rusage[tmp={resources.tmp}] span[ptile={threads}]\" \
    -o \"logs/{rule}.{wildcards}.out\" \
    -e \"logs/{rule}.{wildcards}.err\""
