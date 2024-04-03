# python /net/home/msohail/synergistic_epistasis/scratch2/SFS/src/mac_count_data_byloci.py $IDFILE functype inpath outpath

## This script takes as input
##		input1 (stored in inpath): ALL.dacx.frq.count.indiv.functype.uppercaseaa
##		input2: id_file.txt (a text file with one id name per line)
## This script writes out a file with each individual's name and his/her derived allelic count, one per line. 
## One such file is written for each derived allele count (dac) and functype
## output (stored in outpath): ALL.dacx.frq.count.functype.uppercaseaa.hist

import sys

# read inputs from command line
idfname=sys.argv[1]
inpath = sys.argv[2]
group = sys.argv[3]
level = sys.argv[4]
functype=sys.argv[5]
outpath = sys.argv[6]
status = ""

# create log file 
logf = open(outpath + group + ".dac_hist_byloci.log","a")

# read id file and make hash table
id_dict={}
idf = open(idfname, "r")
ids = idf.readlines()
for i in ids:
    id_dict[i.strip()] = []

# for input derived allele count, check if .frq.count.indiv.functype.uppercaseaa exists. If so parse the file, get 
# 			individual counts and populate the individual dictionary
num_indivs = len(id_dict.keys())
dacs = range(1,2*num_indivs+1)
#dacs = range(1,20)


for dac_i in dacs:
    dac = str(dac_i)
    outf = open(outpath + group + ".ALL.dac" + dac + ".frq.count." + functype +"."+level+ ".uppercaseaa.hist" + status, "w")
    infname = inpath + group + ".dac"+ dac + ".frq.count.indiv." + functype + "."+level+".uppercaseaa" + status
    try:
      # if DAF file found
        #print "SUCCESS"
        inf = open(infname, "r")
        logf.writelines("DAF" + dac + "\t" + functype + "\t" + "PASS" + "\n")
      
        data = inf.readlines()
        #print data
        #counter = 0
        for line in data:
            #counter += 1 
            #scan each line and grab all individuals that have a variant at that line
            indivs = line.strip().split('\t')[7:]
            indivs_names = []
            indivs_counts = []
            #print indivs
            # populate the dictionary
            for indiv in indivs:
              temp = indiv.split(':')
              indivs_names.append(temp[0])
              indivs_counts.append(temp[1])
            #print indivs_names
            #print indivs_counts
            for idii in id_dict.keys():
                if idii in indivs_names:
                    id_dict[idii].append(indivs_counts[indivs_names.index(idii)])
                else:
                    id_dict[idii].append("0")
        inf.close()
        #print id_dict
        #print counter

        # write out final dictionary to outfile
        if len(data)==0:
                for idkey in ids:
                    outf.write("0\n")
        else:	
            for idkey in ids:
                count = id_dict[idkey.strip()]
                #outf.write(idkey.strip())
                for ii in range(len(count)):
                        outf.write(str(count[ii])+" ")
                outf.write("\n")
                id_dict[idkey.strip()] = []

    except:
        #print "FAIL"
	# DAF file not found
          logf.writelines("DAF" + dac + "\t" + functype + "\t" + "FAIL" + "\n")

          for idkey in ids:
              outf.write("0\n")
        
    pass


logf.close()
outf.close()

