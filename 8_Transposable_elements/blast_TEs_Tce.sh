#!/bin/bash

##########################################
#blasting TElibs against Tce TE lib Nosil#
##########################################

module add Blast/ncbi-blast/2.7.1+;

makeblastdb -in Tc_TE_merged_library_renamed_20131009_3_modRM_strict_conv.fa -dbtype nucl -out TceTElibNosil

#blast
for f in *_WickerH.fa; do blastn -task blastn -db TceTElibNosil -query $f -out $f''.blastn -outfmt "6 qseqid sseqid pident qlen slen length mismatch gapopen evalue bitscore" -num_threads$
#for f in *_WickerH.fa; do tblastx -db TceTElibNosil -query $f -out $f''.tblastx -outfmt "6 qseqid sseqid pident qlen slen length mismatch gapopen evalue bitscore" -num_threads 48 -evalu$

module rm Blast/ncbi-blast/2.7.1+;

