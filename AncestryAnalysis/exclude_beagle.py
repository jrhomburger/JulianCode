#! /usr/bin/python

# Script reads in a beaglefile and a file of individuals to exclude
# Outputs a new beaglefile without those individuals

# usage python exclude_beagle.py -b <beaglefile> -e <excludefile> -v vitfile -o <outfile>

import sys
import argparse
import re
parser = argparse.ArgumentParser()
parser.add_argument("-b", help="beagle file name")
parser.add_argument("-e", help="exclude/include individuals file")
parser.add_argument("-include", action="store_true", help="flag, add to include instead of exclude")
parser.add_argument("-v", default="", help="viterbi file name, include to run with a viterbi file")
parser.add_argument("-o", default="outfile", help="outfile_prefix")
parser.add_argument("-ignoreAB", action="store_true", help="tells the program to ignore trailing _A and _B when performing matches on the data")
args=parser.parse_args()

excludefile = open(args.e, "r")
beaglefile = open(args.b, "r")
outfile = open(args.o + ".beagle", "w")

if (args.include):
	print "Including individuals from Beagle File"
else:
	print "Excluding Individuals from Beagle File"

exclude = set()
for elin in excludefile:
	thisind = elin.strip()
	if (args.ignoreAB):
		thisind = re.sub(r'_[AB]$', '', thisind)
	exclude.add(thisind)

beagleheader = beaglefile.readline().strip().split()
if (args.include):
	keep_vec = [0]*len(beagleheader)
else:
	keep_vec = [1]*len(beagleheader)

## Keep header and front of file
keep_vec[0] = 1
keep_vec[1] = 1

if(args.ignoreAB):
	beagleheader = [ re.sub(r'_[AB]$', '' , z) for z in beagleheader]

for k in range(len(beagleheader)):
	if (args.include):
		if (beagleheader[k] in exclude):
			keep_vec[k] = 1
	else:
		if (beagleheader[k] in exclude):
			keep_vec[k] = 0

print "After filtering, " + str(sum(keep_vec)-2) + " of " + str(len(keep_vec)-2) + " haplotypes are included"

headout = [ beagleheader[j] for j in range(len(beagleheader)) if keep_vec[j] == 1]
print "Writing beagle output to " + args.o + ".beagle"
outfile.write("\t".join(headout) + "\n")

for bline in beaglefile:
	bsplit = bline.strip().split()
	bout = [ bsplit[i] for i in range(len(bsplit)) if keep_vec[i] == 1 ]
	outfile.write("\t".join(bout) + "\n")
	
if args.v != "":
	print "Extracting Individuals from associated Viterbi file"
	print "Writing Viterbi output to " + args.o + ".vit"

	vfile = open(args.v, "r")
	vout = open(args.o + ".vit", "w")
	for vline in vfile:
		vsplits = vline.strip().split()
		if (args.ignoreAB):
			thisind = re.sub(r'_[AB]$', '',vsplits[0])
			vsplits[0] = thisind
		else:
			thisind =vsplits[0]
		if (args.include):
			if thisind in exclude:
				 vout.write("\t".join(vsplits) + "\n")
		else:
			if not thisind in exclude:
				 vout.write("\t".join(vsplits) + "\n")
	
	
	
	
	
	
	
	
