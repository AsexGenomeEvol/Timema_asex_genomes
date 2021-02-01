#!/usr/bin/python3.5

import sys

path="cleaned/"
path2="alignment_guidance/final/" 

if len(sys.argv) != 2 :
    sys.exit("Usage: ./last_step.py <prefix: FAM000003>")


main = []
for line in open(path + sys.argv[1] + ".aligned_Hmm10_clean.fasta") :
    if line[0] == ">" :
        main.append(line.rstrip())
        main.append("")
    else :
        main[-1] += line.rstrip()


substitute = []     # only used to get complete header names
for line in open(path2 + sys.argv[1] + ".aligned.fasta") :
    if line[0] == ">" :
        substitute.append(line.rstrip())
        substitute.append("")
    else :
        substitute[-1] += line.rstrip()


if len(main) != len(substitute) :
    sys.exit("Problem with input files (not the same number of sequences)! " + sys.argv[1])


sortie=open(sys.argv[1] + ".final.fasta", "w")
seq_removed=[]
set_aa={'A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V','U','O'}

aln_size=len(main[1])
for i in range(1,len(main),2) :
    if len(main[i]) != aln_size :
        sys.exit(sys.argv[1] + " is not an alignment!")
    
    seq_size = 0
    for j,site in enumerate(main[i]) :
        if site in set_aa :
            seq_size += 1

    if seq_size >= 0.25 * aln_size :
        sortie.write(substitute[i-1] + "\n" + main[i] + "\n")
    else :
        seq_removed.append(main[i-1])
    

sortie.close()
print("\nOutput in " + sys.argv[1] + ".final.fasta!")
print(str(len(seq_removed)) + " sequence(s) were removed !")
if len(seq_removed) > 0 :
    print(str(seq_removed))




