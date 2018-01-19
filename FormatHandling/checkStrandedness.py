### Check strandedness with UCSC ###
## by JRH
## October 24, 2015

## Goal: Read in a bim file and a dbsnp bedfile
## Match rsIDs. Convert position to the position found in the dbsnp bedfile
## Output positions that are not found in the dbSNP file


import argparse
import re

parser=argparse.ArgumentParser()
parser.add_argument("-b", help="Bimfile to convert")
parser.add_argument("-d", help="dbSNP/UCSC file with rsIDs and reference alleles")
parser.add_argument("-r", type=str, help="name of SNPs to flip")
parser.add_argument("-e", type=str, help="name of SNPs that are ambiguous")

args = parser.parse_args()

if (not args.r):
	outfile = args.b + "_reverse"
else:
	outfile = args.r
	
if (not args.e):
	ambifile = args.b + "_ambiguous"
else:
	ambifile = args.e

print "Reading bimfile..."
bimpos = dict()
bimfile = open(args.b)
order = list()
positions_ambi = set()

set_ambi1 = set(["A", "T"])
set_ambi2 = set(["G", "C"])
for b in bimfile:
	bsplit = b.strip().split()
	bimpos[bsplit[1]] = bsplit
	order.append(bsplit[1])
	thisset = set([bsplit[4], bsplit[5]])
	if thisset == set_ambi1 or thisset == set_ambi2:
		positions_ambi.add(bsplit[1])
	
	
dbsnp = open(args.d)
positions_flip = set()
snps_checked = set()
print "Reading UCSC file..."
for d in dbsnp:
	dsplit = d.strip().split()
	if not dsplit[3] in bimpos:
		continue
	snps_checked.add(dsplit[3])
	ref_allele = dsplit[5]
	if not (ref_allele == bimpos[dsplit[3]][4] or ref_allele == bimpos[dsplit[3]][5]):
		positions_flip.add(dsplit[3])
	

### Write the outputs
print "Writing outputs..."
output = open(outfile, "w")

for s in positions_flip:
	output.write(s + "\n")
	
ambi = open(ambifile, "w")

all_snps = set(bimpos.keys())
not_checked = all_snps.difference(snps_checked)

print "Found " + str(len(positions_ambi)) + " ambiguous SNPs"
print "Found " + str(len(not_checked)) + " SNPs not in the UCSC file"

for a in positions_ambi:
	ambi.write(a + "\n")
for b in not_checked:
	ambi.write(b + "\n")


