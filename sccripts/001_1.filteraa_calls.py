# filter_aacalls.py
# python filteraa_calls.py path fname

# This script takes as input:	
#	1) .INFO file (from vcftools --get-INFO AA flag) 
#	2) .frq.count file (from vctools --count2 flag)
# The script outputs:
# 	1) .INFO.uppercaseaa (filtered .INFO file --> contains only those positions with capital letter ancestral allele assignment. Uppercase is when all 3 --(a) human-chimp ancestral
#          sequence, (b) chimp sequence and (c) human-chimp-orangutam ancestral sequences-- all agree )
#	2) .frq.count.uppercaseaa (filtered .frq.count file --> contains only those positions with capital letter ancestral allele assignment)
#	3) .INFO.log (all those positions that have a capital letter ancestral assignment that is diffferent from both extant alleles. These positions are found in the .INFO file but
#          not in the .frq.count file. These positions are not in the final filtered .uppercaseaa files) 

import sys

# read arguements from command line
path=sys.argv[1]
fname=sys.argv[2]

# open input and output files
aaf = open(path + fname + ".INFO", "r")
countf = open(path +  fname + ".frq.count", "r")
aa_outf = open(path + fname + ".INFO.uppercaseaa", "w")
count_outf = open(path + fname + ".frq.count.uppercaseaa", "w")
logfile = open(path + fname + ".INFO.log" , "w")

# read and write out headers
aaheader = aaf.readline()
countheader = countf.readline()
aa_outf.writelines(aaheader)
count_outf.writelines(countheader)

# load all count lines into a dictionary
store = {}
for line in countf.readlines():
    pos = ':'.join((line.split("\t")[0],line.split("\t")[1]))
    store[pos] = line
  
# check each line in .INFO file and write out filtered files with only confident aa calls 
for line in aaf.readlines():
  # grab aa 
    aa = line.strip().split("\t")[4]
    aapos = ':'.join((line.split("\t")[0],line.split("\t")[1]))
  # check aa assignment confidence
    if aa=="A" or aa=="G" or aa=="C" or aa=="T":
      # grab line from .frq.count file
        try:
            countline = store[aapos]
	# write out to filtered files
            aa_outf.writelines(line) 
            count_outf.writelines(countline)
        except:
            errorline = "ERROR: Did not find "  + aapos +  " in .frq.count file\n"
            logfile.writelines(errorline)

aaf.close()
countf.close()
aa_outf.close()
count_outf.close()
