#!/bin/bash
#
# the first argument is the name of species (1_Tdi)
# the second argument is a version indexed reference genome (b3v04)
# the third argument is read group
# (ex: "@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1" )
# note that the output is rather bith bam file, which will be compied to folder data/SP/mapping

#BSUB -L /bin/bash
#BSUB -J "$1"_"$2"_map
#BSUB -q bgee
#BSUB -o ref_is350_to_"$2".log
#BSUB -e ref_is350_to_"$2".err
#BSUB -n 16
#BSUB -M 41943040
#BSUB -R \"rusage[tmp=30000] span[ptile=16]\"

module add UHTS/Aligner/bwa/0.7.17
module add UHTS/Analysis/samtools/1.3

sample="$1"
sp="$2"
reference=data/"$sp"/reference/"$sp"_b3v08.fasta.gz
libs=$(printf -- '%s\n' "${@}" | grep "R1" | cut -f 7,8,9 -d _)
path_to_reads=$(dirname $(printf -- '%s\n' "${@}" | grep "R1" | head -1))

while read index lib; do
    echo Mapping $index $lib
    RG="@RG\tID:"$lib"\tDT:2019-04-23\tLB:lib-"$lib"\tPL:ILLUMINA\tSM:$sample"
    bwa mem -M -t 15 -R \"$RG\" \
        $reference $path_to_reads/"$lib"* | samtools view -bS - > temp.bam
    samtools view -h temp.bam | samblaster -M | samtools view -h -q 20 | samtools sort -@10 -O bam - > "$sample"_"$lib".bam
done < <(printf -- '%s\n' "$libs" | cat -n)

rm temp.bam

ls "$sample"*.bam > "$sample".run_list
samtools merge -@ 16 -b "$sample".run_list ${@:$#}

# rm $(cat "$sample".run_list)