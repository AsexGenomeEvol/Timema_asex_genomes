### TE_landscape_coverage.sh
### want to get coverage of TEs by CpG Kimura distance estimate

## for this we need:
## 1. bam files with all reads aligned (see map_reads.sh)
## 2. alignment files from Repeatmasker (see repeat_masking_local.sh)


# CpG Kimura distances are given in the alignment files. These ests come directly from repeatmasker for versions >4.05 (we used 4.1), and so don't need to correct them with calcDivergenceFromAlign.pl
# I checked by running the align files through calcDivergenceFromAlign.pl - stay the same and same as in the cat file

################################################################
### To get counts for TEs with CpG Kimura distances in the alignment files I need to make a new gff from the align file directly.
### run on Repeatmasker alignment files

for a in algn_files/*.fasta.align; do
echo $a
prefix=`echo $a | sed 's/.fasta.align//' | sed s'/.*\///'`
echo $prefix
python3 TE_align_to_gff.py -a $a -g $prefix -o $prefix
done

### then run HTseq count on mapped read bams

module load gcc
module load htseq/0.11.2

for b in TE_genome_mapping/map/*_aln.sorted.bam; do
echo $b
prefix=`echo $b | sed 's/_aln.sorted.bam//' | sed s'/.*\///'`
echo $prefix
gff_name=`echo $prefix"_b3v08_align.gff" `
echo $gff_name
out_file=`echo $prefix"_htseq_unq_withUNIQID_from_align.counts"`
echo $out_file
htseq-count -f bam -r name -s no -t similarity -i Target --nonunique none $b  $gff_name  > $out_file
done
	

#### add CpG Kimura distance and feature lengths to ests to counts

for c in ./*htseq_unq_withUNIQID_from_align.counts; do
echo $c
prefix=`echo $c | sed 's/.htseq_unq_withUNIQID_from_align.counts//' | sed s'/.*\///'`
echo $prefix
gff_name=`echo $prefix"_b3v08_align.gff" `
echo $gff_name
echo ""
python3 TE_coverage_from_align.py -g $gff_name -c $c -o $prefix
done


### plot results with plot_TE_cov_align.R