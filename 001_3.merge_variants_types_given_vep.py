# -*- coding: utf-8 -*
# merge_id_types.py
# run command: python /net/home/msohail/synergistic_epistasis/src/merge_variants_types_vep.py path fname
#
# 	Input files:
#	1) types = dref0.chr?variant_types.txt
#	chr:pos	mutation_type
#	2) ids = fname.frq.count.indiv.uppercaseaa
#	CHROM   POS     N_ALLELES       N_CHR   {COUNT} indivs
#	
# 	Output files: 
#	1) out =  fname.frq.count.indiv.functype.uppercaseaa
#	CHROM   POS     N_ALLELES       N_CHR   {COUNT} functype indivs
#


import os
import sys

def process(types, ids, out):
	# load everything from ids into a dictionary
	id_dict = {}
	fids = open(ids, 'r')
	for line in fids.readlines():
		split = line.split("\t")
		posid = ":".join(split[0:2])
		#print posid
		#print line
		id_dict[posid] = line
	#	print split[0]
	#	print split[1]
	fids.close()


	# read in types file and write out id_types file
	ftypes = open(types, 'r')
	fout = open(out, 'w')
	for line in ftypes.readlines():
		split = line.split("\t")
               # print(split)
		pos = split[0].strip()
		mtype = split[1].strip()
	#	print pos
	#	print mtype
		if pos in id_dict:
		  curr_line = id_dict[pos].strip()
		  curr_split = curr_line.split("\t")
		  indivs = "\t".join(curr_split[6:])
		  printline = "\t".join(curr_split[0:6])
		  outline = printline + "\t" + mtype + "\t" + indivs + "\n"
		  fout.writelines(outline)
		else:
		  #print "ERROR " , pos, " not found in " + ids 
		  pass
	ftypes.close()
	fout.close()


## ENTER INPUT AND OUTPUT FILE PATHS HERE

# read in and set names and prefix to append at end of file names
path = sys.argv[1]
group = sys.argv[2]
types =  path + group + ".S.INFO.txt" # types: vep -> S
ids = path + group+".frq.count.indiv.uppercaseaa"
out =  path + group+".frq.count.indiv.functype.S.uppercaseaa"


process(types, ids, out)

