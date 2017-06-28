module add UHTS/Analysis/samtools/1.3

cp /archive/dee/schwander/ptranvan/mites/assembly/spades/without_conta/sm2/round3/sspace/111/Gapcloser/bwa/Sm2.trim.{180,350,550}.pair12.sam .

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

sort_and_assign_RG Sm2 180 r1 &
sort_and_assign_RG Sm2 350 r2 &
sort_and_assign_RG Sm2 550 r3 &

wait

rm Sm2.trim.{180,350,550}.pair12.sam
samtools merge ref_to_2_Sm_b1v01.bam Sm2.trim.{180,350,550}.pair12.sort.RG.sam
samtools index ref_to_2_Sm_b1v01.bam

rm Sm6.trim.{180,350,550}.pair12.sort.RG.sam

atlas task=estimateTheta \
	bam=ref_to_2_Sm_b1v01.bam \
	window=10000 \
	suppressWarnings verbose \
	1> ref_to_2_Sm_b1v01_theta.log
    Sm2.trim.{180,350,550}.pair12.sam
