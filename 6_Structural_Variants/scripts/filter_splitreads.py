#!/usr/bin/env python3
import pysam
import argparse
import sys
import time

def parseArgs():
    parser = argparse.ArgumentParser(description="Filter read pairs that don't map to the same scaffold.")
    parser.add_argument("inBAMfile", help="The BAM file to be filtered")
    parser.add_argument("outBAMfile", help="The BAM file to be generated")
    return parser

def main(args):
    args = parseArgs().parse_args(args)

    start = time.time()
    # input_file = 'data/3_Tms/mapping/Tms_01_to_b3v08.bam'
    in_bamfile = pysam.AlignmentFile(args.inBAMfile, "rb")

    # If I would like to adjust header in here:
    # header = samfile.header.as_dict()
    # # consider adding something to PG, so it's clear that it was not simply produced by bwa
    #
    # for key in header['HD']:
    #      if header[key] < 1000:
    #           remove it from header
    # filtered_bam = pysam.AlignmentFile("data/3_Tms/mapping/Tms_01_to_b3v08_same_scf_reads.bam", "wb", header = header)
    # filtered_bam = pysam.AlignmentFile("data/3_Tms/mapping/Tms_01_to_b3v08_same_scf_reads.bam", "wb", template=samfile)
    out_bamfile = pysam.AlignmentFile(args.outBAMfile, "wb", template=in_bamfile)

    parsed = 0
    kept = 0
    for read in in_bamfile.fetch():
        # read.reference_name
        parsed += 1
        if read.reference_name == read.next_reference_name:
            kept += 1
            out_bamfile.write(read)

    end = time.time()

    sys.stderr.write("Processed " + str(parsed) + " reads.\n")
    sys.stderr.write("Kept " + str(kept) + " of them (" + str(round((kept / parsed) * 100, 2)) + "%).\n")
    sys.stderr.write("   in " + str(end - start) + " seconds\n")

if __name__ == "__main__":
    args = None
    if len(sys.argv) == 1:
        args = ["--help"]
    main(args)