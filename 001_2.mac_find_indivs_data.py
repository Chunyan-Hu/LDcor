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
filter = sys.argv[3]
posf = open(path + fname + ".frq.count" + "." + filter, "r")
gtf = open(path + fname + "." + filter + ".GT.FORMAT", "r")
infof = open(path + fname + ".INFO" + "." + filter, "r")
outf = open(path + fname + ".frq.count.indiv" + "." + filter, "w")

# grab headers and sample ids
headout = posf.readline().strip() + "\tindivs\n"
infoheader = infof.readline()
header = gtf.readline()
ids = header.strip().split("\t")[2:]
outf.writelines(headout)


def is_derived(ref, alt, aa):
	""" inputs: ref = reference allele, alt = alternative allele, aa = ancestral allele
		determine if the derived allele is the reference or the alternative
	 	output: return 'ref' or 'alt' respectively """
	if aa.lower() == ref.lower():
		derived = "alt"
	else:
		derived = "ref"
	return derived
	
def count_derived(check, gt):
	""" inputs: check = 'ref' or 'alt' ; gt = '0|1' or '1|0' etc
		count the number of derived alleles in the genotype
		output: return count """
	# if derived allele is alternative
	if check == "alt" :
		if currgt == "0|1" or currgt == "1|0" or currgt == "0/1" or currgt == "1/0":
		  count = 1
		elif currgt == "1|1" or currgt == "1/1":
		  count = 2
		else: 
		  count = False
	# if derived allele is reference
	else:
		if currgt == "0|1" or currgt == "1|0" or currgt == "0/1" or currgt == "1/0":
		  count = 1
		elif currgt == "0|0" or currgt == "0/0":
		  count = 2
		else: 
		  count = False
	return count

	
# main body	
for line in gtf.readlines():
    gts = line.strip().split("\t")[2:]
    outline = posf.readline().strip()
  
  # get info, determine ancestral/derived
    infoline = infof.readline().strip()
  #print infoline
    parts = infoline.strip().split("\t")
    ref = parts[2]
    alt = parts[3]
    aa = parts[4]
    check = is_derived(ref, alt, aa)
  
  # calculate and write out genotype count for each individual
    for i in range(len(gts)):
        currid = ids[i]
        currgt = gts[i]
        count = count_derived(check, currgt)
      
        if count:
            outline = outline + "\t" + currid + ":"+ str(count)
    outf.writelines(outline+"\n")   
  
outf.close()

posf.close()
gtf.close()
infof.close()

