# python /net/home/msohail/synergistic_epistasis/src/mac_count_data.py $IDFILE functype inpath outpath
# python /net/home/msohail/scratch2/SFS/synergistic_epistasis/src/mac_count_data.py /net/home/msohail/scratch2/SFS/synergistic_epistasis/provinces/twinsuk.o possibly /net/home/msohail/scratch2/SFS/synergistic_epistasis/commonvariants/twinsuk/data/ /net/home/msohail/scratch2/SFS/synergistic_epistasis/commonvariants/twinsuk/histfiles/

## This script takes as input
##		input1 (stored in inpath): ALL.dacx.frq.count.indiv.functype.uppercaseaa
##		input2: id_file.txt (a text file with one id name per line)
## This script writes out a file with each individual's name and his derived allelic count, one per line. 
## One such file is written for each derived allele count (dac) and functype
## output (stored in outpath): ALL.dacx.frq.count.functype.uppercaseaa.hist

import sys

# read inputs from command line
idfname=sys.argv[1]
inpath =sys.argv[2]
group= sys.argv[3]
level=sys.argv[4]
functype=sys.argv[5]
outpath = sys.argv[6]
status = ""

# create log file 
logf = open(outpath +group+ ".dac_hist.log","a")

# read id file and make hash table
id_dict={}
idf = open(idfname, "r")
ids = idf.readlines()
for i in ids:
    id_dict[i.strip()] = 0

# for input derived allele count, check if .frq.count.indiv.functype.uppercaseaa exists. If so parse the file, get 
# 			individual counts and populate the individual dictionary
num_indivs = len(id_dict.keys())
dacs = range(1,2*num_indivs+1)
#dacs = range(1,23)

for dac_i in dacs:
    #print(dac_i)
    dac = str(dac_i)
	#outf = open(outpath + "/" + dac + ".frq.count.uppercaseaa.hist"+ status, "w")
	#infname = inpath + "/" + dac + ".frq.count.indiv.uppercaseaa" + status
    outf = open(outpath + group+".ALL.dac" + dac + ".frq.count." + functype +"."+level+ ".uppercaseaa.hist"+ status, "w")
    infname = inpath + group + ".dac"+ dac + ".frq.count.indiv." + functype + "."+level+ ".uppercaseaa" + status
try:
	  # if DAF file found
    inf = open(infname, "r")
    logf.writelines("DAF" + dac + "\t" + functype + "\t" + "PASS" + "\n")
	  
    data = inf.readlines()
	  #print data
    for line in data:
		#scan each line and grab all individuals that have a variant at that line
        indivs = line.strip().split('\t')[7:]
		#print indivs
		# populate the dictionary
    for indiv in indivs:
        temp = indiv.split(':')
        id_dict[temp[0]] += int(temp[1])
        inf.close()
	  
except:
	  # DAF file not found
       logf.writelines("DAF" + dac + "\t" + functype + "\t" + "FAIL" + "\n")
pass

	# write out final dictionary to outfile
for idkey in ids:
    count = id_dict[idkey.strip()]
    outf.writelines(idkey.strip() + "\t" + str(count) + "\n")
    id_dict[idkey.strip()] = 0

logf.close()
outf.close()
