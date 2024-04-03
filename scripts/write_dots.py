## python ~/scratch2/SFS/synergistic_epistasis/src/write_dots.py 11.INFO.uppercaseaa

import sys
fname = sys.argv[1]
f = open(fname, "r")
outf = open(fname+ ".vepinput", "w")

header = next(f)
tmp = header.split()
outheader = "\t".join(tmp[0:2]) + "\t" + "ID" + "\t" + "\t".join(tmp[2:]) + "\n"
outheader2 = "#" + outheader
outf.writelines(outheader2)
for line in f:
	tmp = line.split()
	outline = "\t".join(tmp[0:2]) + "\t.\t" + "\t".join(tmp[2:]) + "\n"
	outf.writelines(outline)

f.close()
outf.close()
