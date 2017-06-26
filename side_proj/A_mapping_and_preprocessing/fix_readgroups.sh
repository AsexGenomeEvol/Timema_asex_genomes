#!/bin/bash


# echo -e "@RG\tOs5.trim.180.pair12.sort
# @RG\tOs5.trim.350.pair12.sort
# @RG\tOs5.trim.550.pair12.sort" > rg.txt
# samtools view -H Os5.trim.180.pair12.sam >> rg.txt
# samtools merge -rh rg.txt ref_to_1_Os_b1v01.bam Os5.trim.180.pair12.sort.sam Os5.trim.350.pair12.sort.sam Os5.trim.550.pair12.sort.sam
# samtools index ref_to_1_Os_b1v01.bam
# atlas task=estimateTheta \
# 	bam=ref_to_1_Os_b1v01.bam \
# 	window=10000 \
# 	suppressWarnings verbose \
# 	1> ref_to_1_Os_b1v01_theta.log
#
#
# samtools view -H ref_to_1_Os_b1v01.bam > header.sam
# cat header.sam | awk '{if($2 == "ID:bwa"){ printf("@RG\tZ:Os5.trim.180.pair12.sort\n@RG\tZ:Os5.trim.350.pair12.sort\n@RG\tZ:Os5.trim.550.pair12.sort\n%s\n",$0) } else {print $0}}' > header_corrected.sam
# samtools reheader header_corrected.sam ref_to_1_Os_b1v01.bam
#
#
# samtools view -H ref_to_1_On_b1v01.bam > header.sam
#
# samtools reheader header_corrected.sam ref_to_1_On_b1v01.bam

samtools sort Os5.trim.180.pair12.sort.sam Os5.trim.180.pair12.sam
samtools sort Os5.trim.350.pair12.sort.sam Os5.trim.350.pair12.sam
samtools sort Os5.trim.550.pair12.sort.sam Os5.trim.550.pair12.sam

samtools sort - Os5.trim.550.pair12.sam |

function sort_and_assign_RG {
    # $1 sp
    # $2 IS
    # $3 RG ID
    SP=$1
    IS=$2
    RG="ID:"$3"\tSM:ref\tLB:is"$IS
    IFILE=$SP.trim.$IS.pair12.sam
    OFILE=$SP.trim.$IS.pair12.sort.RG.sam
    HEADER=$SP.$IS.header.sam

    echo -e "@HD\tVN:1.3\tSO:coordinate" > $HEADER
    samtools view -H $IFILE | awk -v RG="$RG" '{if($2 == "ID:bwa"){ printf("@RG\t%s\n%s\n", RG, $0) } else {print $0}}' >> $HEADER
    samtools sort -- $IFILE | samtools view - | awk -v RGID="$3" '{ printf "%s\tRG:Z:%s\n",$0, RGID; }' | cat $HEADER - > $OFILE
    rm $HEADER
}

function assign_RG {
    # $1 sp
    # $2 IS
    # $3 RG ID
    SP=$1
    IS=$2
    RG="ID:"$3"\tSM:ref\tLB:is"$IS
    IFILE=$SP.trim.$IS.pair12.sort.sam
    OFILE=$SP.trim.$IS.pair12.sort.RG.sam
    HEADER=$SP.$IS.header.sam

    samtools view -H $IFILE | awk -v RG="$RG" '{if($2 == "ID:bwa"){ printf("@RG\t%s\n%s\n", RG, $0) } else {print $0}}' >> $HEADER
    samtools view $IFILE | awk -v RGID="$3" '{ printf "%s\tRG:Z:%s\n",$0, RGID; }' | cat $HEADER - > $OFILE
    rm $HEADER
}

assign_RG Os5 180 r1 &
assign_RG Os5 350 r2 &
assign_RG Os5 550 r3 &

samtools merge ref_to_1_Os_b1v01.bam $SP.trim.{180,350,550}.pair12.sort.RG.sam
samtools index ref_to_1_Os_b1v01.bam
atlas task=estimateTheta \
	bam=ref_to_1_Os_b1v01.bam \
	window=10000 \
	suppressWarnings verbose \
	1> ref_to_1_Os_b1v01_theta.log

assign_RG On6 180 r1 &
assign_RG On6 350 r2 &
assign_RG On6 550 r3 &
samtools merge ref_to_1_On_b1v01.bam On6.trim.{180,350,550}.pair12.sort.RG.sam
samtools index ref_to_1_On_b1v01.bam

atlas task=estimateTheta \
	bam=ref_to_1_On_b1v01.bam \
	window=10000 \
	suppressWarnings verbose \
	1> ref_to_1_On_b1v01_theta.log

git@github.com:AsexGenomeEvol/timema_assembly.git
