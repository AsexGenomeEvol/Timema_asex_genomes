
```
VARIANTS=test/test_diploidSV.vcf.gz
REF=test/Tms_sample_ref.fasta
TESTOUT=test/paragraph_compatible.vcf
SCRIPT=scripts/convertManta/convertManta2Paragraph_compatible_vcf.py
python3 $SCRIPT $REF $VARIANTS > $TESTOUT
```

This will be used for testing of genotyping (should sub-select only the first 10 chromosomes though)


### Generic genomics

scripts I with wider possible use. They are developed in separated repository that is linked here using submodule.
