# python /net/home/msohail/scratch2/SFS/synergistic_epistasis/src/mac_find_indivs_data.py /net/data/msohail/commonvariants/ALS_illumina/data/ 22 uppercaseaa

## This script takes the output of vcftools --counts2, --extract-FORMAT-info GT and --INFO aa flags as input
##	input: fname.frq.count.filter
##		   fname.filter.GT.FORMAT
##		   fname.INFO.filter
## This script outputs a new file which appends individual derived allele counts at each position
## output: fname.frq.count.indiv.filter

import sys

# open input and output files
path=sys.argv[1]
fname=sys.argv[2]
#filter = sys.argv[3]
gtf = open(path + "/" +  fname + ".uppercaseaa.GT.FORMAT", "r")
infof = open(path + "/" + fname + ".INFO.uppercaseaa", "r")
outf = open(path + "/" +  fname + ".alt.count.indiv", "w")

# grab headers and sample ids
infoheader = infof.readline()
headout = str(infoheader.strip().split("\t")[:4])+ "\tindivs\n"
header = gtf.readline()
ids = header.strip().split("\t")[2:]
outf.writelines(headout)


def count_alt(gt):
    if currgt == "0|1" or currgt == "1|0" or currgt == "0/1" or currgt == "1/0":
        count = 1
    elif currgt == "1|1" or currgt == "1/1":
        count = 2
    elif currgt == "0|0" or currgt == "0/0":
        count = 0
    else: 
        count = 9
    return count

	
# main body	
for line in gtf.readlines():
    #print(line.strip().split("\t")[1:2])
    gts = line.strip().split("\t")[2:]
    #print(gts[1])
  # get info, determine ancestral/derived
    infoline = infof.readline().strip()
    #print(infoline)
    parts = infoline.strip().split("\t")
    #print(len(parts))
    ref = parts[2]
    alt = parts[3]
    outline = infoline
  # calculate and write out genotype count for each individual
    for i in range(len(gts)):
        #currid = ids[i]
        currgt = gts[i]
        count = count_alt(currgt)
        #if count:
        outline = outline + "\t" + str(count)
    outf.writelines(outline+"\n")   
  
outf.close()

gtf.close()
infof.close()

