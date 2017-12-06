#!/usr/bin/python3.5

from __future__ import print_function
import sys
import os
import math


# Script to manually filter positions using JEXL-like expressions !
# (as there are issues with GATK 'VariantFiltration', hard filtering in the FORMAT field is not applied at every position...) 
# Script works on : 
#           -> SNP/indel vcf
#           -> all position vcf (as well)

# Format of the option file :
# example :
#-----------
#INFO:QD<5.0:lowQD
#INFO:DP>100:highDP
#FORMAT:DP<10:lowDP
#FORMAT:GQ<30:lowGQ
#---------
# type of the field ('FORMAT' or 'INFO') should be given first, 
#      then the JEXL expression (with no space!),
#      then FILTER expression to write if position fails the test.
# use ':' as separator !
# threshold values can be given as int or float, they will converted to float !



# Check input :
if len(sys.argv) != 4 :
    if len(sys.argv) != 5 :
        print("Usage: ./3_variant_filtration_r1.py  <vcf_file>  <headers_file>  <criteria_file>  (<NO_AMB>)")
        print("\n       -> 'NO_AMB' tag (at the end) : optional parameter to filter out ambiguous (ie indels & multiallelic) sites")
        sys.exit("                                      use it (if wanted) on SNP or allSites vcf, but not on INDEL vcf...")
    if len(sys.argv) == 5 and sys.argv[4] != "NO_AMB" :
        print("Usage: ./3_variant_filtration_r1.py  <vcf_file>  <headers_file>  <criteria_file>  (<NO_AMB>)")
        print("\n       -> 'NO_AMB' tag (at the end) : optional parameter to filter out ambiguous (ie indels & multiallelic) sites")
        sys.exit("                                      use it (if wanted) on SNP or allSites vcf, but not on INDEL vcf...")
    


# Read input :
vcf_file =   sys.argv[1]    # it will be read line by line
headers  = [ x.rstrip() for x in open(sys.argv[2]) if x[0]=="#" ]
filters  = [ x.rstrip().split(":") for x in open(sys.argv[3]) ]

for line in filters :
    if len(line) != 3 :
        sys.exit("\n\nERROR: filtering file badly formatted !")

accepted_FORMAT_param = {'QD','FS','SOR','MQ','MQRankSum','ReadPosRankSum','DP'}
accepted_INFO_param   = {'DP','GQ'}

# Get number of ind :
if headers[-1].split()[0] != "#CHROM" :
    sys.exit("\n\nERROR: vcf badly formatted, issue with the last header line : '#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT ...'")
nb_inds = len(headers[-1].split()) - 9

# Filter out ambiguous sites :
filter_amb = False
if len(sys.argv) == 5 :
    filter_amb = True

# Print filters (so they can also appear in a log) :
print("\n\nHard filtering parameters, thresholds and associated tags :\n(even if FORMAT tags will be merged into a single annotation)\n")
for line in filters :
    if line[0] == 'INFO' :
        line = line[0]+":    "+line[1]+"  (tag:'"+line[2]+"')"
    elif line[0] == 'FORMAT' :
        line = line[0]+":  "+line[1]+"  (tag:'"+line[2]+"')"
    print(line)



#----------
# FUNCTIONS :

def compare(a,operand,b) :
    # do arithmetic comparisons :
    if operand == ">=" :
        if a >= b :
            return True
        else :
            return False
    if operand == "<=" :
        if a <= b :
            return True
        else :
            return False
    if operand == "=" :
        if a == b :
            return True
        else :
            return False
    if operand == ">" :
        if a > b :
            return True
        else :
            return False
    if operand == "<" :
        if a < b :
            return True
        else :
            return False



def INFOfiltering(pos,filters) :

    jexl = filters[0]
    expr = filters[1]
    
    operands = [">=","<=","=",">","<"]

    # decompose jexl expression :
    operand = ""
    for op in operands :
        if op in jexl :
            operand = op
            jexl = jexl.replace(operand,":").split(":")
            param = jexl[0]
            thres = float(jexl[1])
            break
    if operand == "" :
        sys.exit("\n\nERROR: bad format for jexl expression in line : '"+str(jexl)+"' of filtering file !")

    # find parameter value in INFO field :
    INFO = pos[7].split(";")
    valid = False
    for p in INFO :
        if p.split("=")[0] == param :
            value = float(p.split("=")[1])     # float("NaN") -> nan (not a number), 'MQ' sometimes has 'NaN' value..
            if compare(value,operand,thres) or math.isnan(value) :   # position fails the test
                if pos[6] == "." :
                    pos[6] = expr
                else :
                    pos[6] += ";"+expr
            else :
                valid = True
    return pos[6],valid



def FORMATfiltering(pos,filters) :

    jexl = filters[0]
    expr = filters[1]
    
    operands = [">=","<=","=",">","<"]

    # decompose jexl expression :
    operand = ""
    for op in operands :
        if op in jexl :
            operand = op
            jexl = jexl.replace(operand,":").split(":")
            param = jexl[0]
            thres = float(jexl[1])
            break
    if operand == "" :
        sys.exit("\n\nERROR: bad format for jexl expression in line : '"+str(jexl)+"' of filtering file !")

    # find parameter position in FORMAT field :
    FORMAT = pos[8].split(":")
    FORMAT_len = len(FORMAT)
    found = False
    param_pos=-1

    # look for the parameter position in FORMAT field :
    for i,p in enumerate(FORMAT) :
        if p == param :
            found = True
            param_pos = i
            break
    if param_pos == -1 :     # not found
        warn = "not_found"
        return warn,found

    # check the correct number of ind :
    if len(pos) - 9 != nb_inds :
        # warn = "incorrect_ind_nb"
        # return warn,found
        sys.exit("\n\nERROR: missing sample at line :\n"+"   ".join(pos))

    # search for the parameter in each ind :
    result_per_ind = []
    
    for ind in pos[9:] :
        ind = ind.split(":")

        if len(ind) != FORMAT_len :     # check that all inds have the correct number of sub-fields
            # warn = "incorrect_subfield_nb"
            # return warn,found
            sys.exit("\n\nERROR: incorrect number of FORMAT subfields at line :\n"+"   ".join(pos))  
            # should not happen anymore as FORMAT subfields were completed !

        value = ind[param_pos]
        if value == '.' or value == '.,.' or value == '.,.,.' or value == '.|.' :
            result_per_ind.append('.')
            continue
        value = float(value)
        if compare(value,operand,thres) == True :  # ind fails this test
            result_per_ind.append(expr)
        else :
            result_per_ind.append('pass')

    if len(result_per_ind) != nb_inds :
        sys.exit("\n\nERROR: see DEBUG_001 in script !")

    return result_per_ind,found



def FORMATcomplete(pos) :

    FORMAT = pos[8].split(":")
    FORMAT_len = len(FORMAT)      # no 'FT' subfield at this step
    
    set1 = {'DP','GQ','RGQ','PS','GP','PQ','PID'} # -> '.'
    set2 = {'AD','HQ'} # -> '.,.'
    set3 = {'GL','PL'} # -> '.,.,.'
    set4 = {'PGT'}     # -> '|'

    for i,ind in enumerate(pos[9:]) :
        ind = ind.split(":")
        
        if len(ind) > FORMAT_len :
            sys.exit("\n\nERROR: too many subfields in sample "+str(i+1)+" :\n"+"   ".join(pos))
        elif len(ind) < FORMAT_len :  # if sample is missing subfields
            ind_new = ind[:]
            for j in range(len(ind),FORMAT_len) :
                if FORMAT[j] in set1 :
                    ind_new.append(".")
                elif FORMAT[j] in set2 :
                    ind_new.append(".,.")
                elif FORMAT[j] in set3 :
                    ind_new.append(".,.,.")
                elif FORMAT[j] in set4 :
                    ind_new.append(".|.")
                else :
                    ind_new.append(".")
                    print("WARNING: subfield "+FORMAT[j]+" was completed with '.'!")

            pos = pos[:9+i] + [":".join(ind_new)] + pos[9+i+1:]
    return pos





#----------
# MAIN LOOP :

output_name_content = sys.argv[1].split("/")[-1].replace(".vcf",".filtered.vcf")     # remove path
output_content = open(output_name_content, "w")   # will not contain headers !

# loop over vcf positions :
i=-1
with open(vcf_file) as vcf :
    for pos in vcf :
    
        i += 1
        if pos[0] == "#" :    # in case headers were not removed (they are given in a separate file)
            continue
        else :
            pos = pos.rstrip().split()
    
        if i+1 % 1000000 == 0 :       # not working for unknown reason..
            print("   -> "+str(i+1)+" positions done !")
            sys.stdout.flush()
    
        validatedF = False      # to be validated, a position should fail no filter, 
                                # and have the information for at least one of them.
    
        FORMAT_results = []     # will contain for each ind, the results per parameter.
        for ind in range(0,nb_inds) :
            FORMAT_results.append([])
    

        # complete missing subfields in samples : (if there are some)
        pos = FORMATcomplete(pos)


        #-------------------------------
        # loop over filtering parameters :
        for j,param in enumerate(filters) :
    
            # INFO parameter :
            if param[0] == 'INFO' :
                pos[6], valid = INFOfiltering(pos, param[1:])
                if valid == True :
                    validatedF = True
    
            # FORMAT parameter :
            elif param[0] == 'FORMAT' :
                output, found = FORMATfiltering(pos, param[1:])
                if output == "not_found" : # or output == "incorrect_ind_nb" or output == "incorrect_subfield_nb" : # last two cases produce ERROR !
                    output = nb_inds * ['.']
                
                for ind in range(0,nb_inds) :
                    FORMAT_results[ind].append(output[ind])
    
            else :
                sys.exit("\n\nERROR: only 'FORMAT' & 'INFO' fields can be used (not '"+str(thres[0])+"'), modify your input filtering file !")
            #---------------------------------------------
    
    
        # Has the position PASS the INFO tests :
        amb_site = False
        if filter_amb == True :
            if len(pos[3]) > 1 or len(pos[4]) > 1 or pos[4] == '*' :
                amb_site = True

        if amb_site == False :
            if validatedF == True and pos[6] == '.' :
                pos[6] = 'PASS'
            if validatedF == False and pos[6] == '.' :      # ie, no field was found for any of the parameter tested...
                pos[6] = 'x' 
        else :
            pos[6] = 'ambiguous'    # we remove ambiguous (indel/multiallelic) sites
    
    
        # Has the position PASS the FORMAT tests :
        
        FORMAT_order = sorted(pos[8].split(":")[1:])    # small check on subfields order
        if pos[8].split(":")[1:] != FORMAT_order :             
            sys.exit("\n\nERROR: FORMAT subfields are not in alphabetical order in ('GT' always first) :\n"+"   ".join(pos))
        
        pos[8] = pos[8]+":FT"
        FORMAT_subfields = sorted(pos[8].split(":"))
        if "GT" in FORMAT_subfields :
            FORMAT_subfields.remove("GT")
            FORMAT_subfields = ["GT"] + FORMAT_subfields    # 'GT' field should always appear first (if present).
        pos[8] = ":".join(FORMAT_subfields)
        for k,subfield in enumerate(FORMAT_subfields) :
            if subfield == "FT" :
                FT_pos = k
                break
   
        # loop on samples :
        for ind in range(0,nb_inds) :
            
            # output formatting :
            sample_begin = ":".join(pos[9+ind].split(":")[:FT_pos]) + ":"  # 'FT' will never be first field (-> 'GT') !
            if len(pos[9+ind].split(":")[FT_pos:]) > 0 :
                sample_end = ":" + ":".join(pos[9+ind].split(":")[FT_pos:])
            else :
                sample_end = ""
       

            if pos[6] != 'PASS' :     # sample can pass only if the SNP was validated at the population level.
                pos[9+ind] = sample_begin+"x"+sample_end
    
            # sample has succeed at least one test and failed none :
            elif set(FORMAT_results[ind]) == {'pass'} or set(FORMAT_results[ind]) == {'pass','.'} :
                pos[9+ind] = sample_begin+"PASS"+sample_end
    
            # sample could not be tested for any of the parameters ('x' tag) :
            # (either because of missing field in FORMAT or value in sample not informed, ie equals '.')
            elif set(FORMAT_results[ind]) == {'.'} :
                pos[9+ind] = sample_begin+"x"+sample_end

            # sample failed at least one test ('f' tag) :
            else :
                # we do not write complete list of failed tags anymore as it significantly increase file size !
                # tag = set(FORMAT_results[ind])
                # tag.discard('pass')
                # tag.discard('.')
                # tag = ";".join(sorted(list(tag)))
                # pos[9+ind] = ":".join(pos[9+ind].split(":")[:FT_pos])+":"+tag+":"+":".join(pos[9+ind].split(":")[FT_pos:])
                pos[9+ind] = sample_begin+"f"+sample_end
    
        # write output :
        output_content.write("\t".join(pos)+"\n")

# end of main loop
#-----------------




# Write header file (if necessary) :

name_headers        = sys.argv[2]
output_name_headers = name_headers+".filtered"

found_head = False
if os.path.isfile(output_name_headers) :    # we will write only if it has not been created !
    found_head = True
else :
    output_headers = open(output_name_headers, "w")

    # complete '##FILTER' lines :
    new_head_FILTER = []
    line_FILTER = -1
    for i,line in enumerate(headers) :
        if line == '##FILTER=<ID=LowQual,Description="Low quality">' :
            line_FILTER = i
            for param in filters :
                if ">=" in param[1] or "<=" in param[1] :
                    param[1] = param[1].replace("<="," <= ").replace(">="," >= ")
                else :
                    param[1] = param[1].replace("<"," < ").replace(">"," > ").replace("="," = ")
                new_head_FILTER.append('##FILTER=<ID='+param[2]+',Description="'+param[1]+'">')
            break
    if line_FILTER == -1 :
        sys.exit("\n\nERROR: output file headers could not be updated, '##FILTER=<ID=LowQual,Description=\"Low quality\">' line missing !")
    else :
        headers = headers[:line_FILTER+1] + new_head_FILTER + headers[line_FILTER+1:]
     
    # complete '##FORMAT' line :
    new_head_FORMAT = []
    line_FORMAT = -1
    for i,line in enumerate(headers) :
        if line == '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">' :
            line_FORMAT = i
            new_head_FORMAT.append('##FORMAT=<ID=FT,Number=.,Type=String,Description="Genotype-level filter">')
            break
    if line_FORMAT == -1 :
        sys.exit("\n\nERROR: output file headers could not be updated, '##FORMAT=<ID=DP,Number=1,Type=Integer,[...]' line missing !")
    else :
        headers = headers[:line_FORMAT+1] + new_head_FORMAT + headers[line_FORMAT+1:]

    # write header output :
    for i,line in enumerate(headers) :
        output_headers.write(line+"\n")


# print info in console :

print("\nList and meaning of tags :")
print("FORMAT field :")
print("               * 'PASS'   (position has passed at least one test and failed none)")
print("               *  list of failed tags (ex: 'lowQD;lowMQ')")
print("               *  'x'     (no test could be performed because of missing fields, not expected)")
print("               * 'ambiguous'  (if NO_AMB option was set, meaning position is indel or multiallelic)")
print("INFO field :")
print("               * 'PASS' (sample has passed at least one test and failed none)")
print("               * 'f'    (sample has failed at least one of the test)")
print("               * 'x'    (the position failed the FORMAT filter, then all samples will get this tag)")
print("               *        (or no test could be performed : missing field(s) in INFO or value(s) undefined for specific sample)")

print("\nOutput for '"+vcf_file+"' in '"+output_name_content+"' (no header in this output file) !")
if found_head == False :
    print("Output for '"+name_headers+"' in '"+output_name_headers+"' !")
else :
    print("\nWARNING: '"+output_name_headers+"' already exists (this is expected with parallelisation) !")



