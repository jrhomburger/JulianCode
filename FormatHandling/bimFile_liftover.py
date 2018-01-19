### Update hg18 to hg19 ###
## by JRH
## October 24, 2015

## Goal: Read in a bim file and a dbsnp bedfile
## Match rsIDs. Convert position to the position found in the dbsnp bedfile
## Output positions that are not found in the dbSNP file


import argparse
import re

parser=argparse.ArgumentParser()
parser.add_argument("-b", help="Bimfile to convert")
parser.add_argument("-d", help="dbSNP file with correct positions")
parser.add_argument("-o", type=str, help="outputfile name, default is bimfile_updated")
parser.add_argument("-e", type=str, help="name of snps that couldn't be mapped")

args = parser.parse_args()

if (not args.o):
	outfile = args.b + "_updated"
else:
	outfile = args.o

print "Reading bimfile"
bimpos = dict()
bimfile = open(args.b)
order = list()
for b in bimfile:
	bsplit = b.strip().split()
	bimpos[bsplit[1]] = bsplit
	order.append(bsplit[1])
	
dbsnp = open(args.d)
positions = dict()
outdict = dict()
print "Reading dbSNP file"
for d in dbsnp:
	if (d.startswith("track")):
		continue
	dsplit = d.strip().split()
	if len(dsplit) < 4:
		continue
	if not dsplit[3] in bimpos:
		continue
		
	if (len(dsplit[0]) == 4):
		chr = dsplit[0][-1]
	elif (len(dsplit[0]) == 5):
		chr = dsplit[0][-2:]
	pos = dsplit[2]
	outdict[dsplit[3]] = [chr, dsplit[3], "0", pos, bimpos[dsplit[3]][4], bimpos[dsplit[3]][5]]
	
	
print "Writing output"
out = open(outfile, "w")
if (args.e):
	excludefile=open(args.e, "w")

for o in order:
	if o in outdict:
		out.write("\t".join(map(str, outdict[o])) + "\n")
	else:
		out.write("\t".join(map(str, bimpos[o])) + "\n")
		excludefile.write(bimpos[o][1] + "\n")

