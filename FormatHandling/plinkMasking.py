### Need to generate sites to zero out in masking data based upon the bed file

import argparse
import sys
import os.path
import re

parser = argparse.ArgumentParser()

parser.add_argument("-bim", help="bimfile of PLINk file to mask")
parser.add_argument("-inds", help="list of individuals to mask, with a FamID then a space/tab and IndID")
parser.add_argument("-onlyCluster", help="if turned on, only writes out cluster mask files (doesn't run through PLINK")
parser.add_argument("-diffNames", help="if turned on, looks for the bedfile name for each individual as the third column of the indfile, if not bedfile name will be sampID.bed")
parser.add_argument("-keep", help="name of region in the bedfile to keep, or name of file with names of regions to keep")
parser.add_argument("-cstart", default=1, type=int, help="number of the cluster to start on")

args=parser.parse_args()

keep = set()
## read in the masking strings
if (os.path.isfile(args.keep)):

	masksfile = open(args.keep, "r")
	for m in masksfile:
		keep.add(m.strip())
else: 
	keep.add(args.keep)

## First, read in the bimfile and all of the positions

print "Reading bimfile..."
all_pos = dict() # chromosome indexed dictionary of positions
bimfile = open(args.bim, "r")

for b in bimfile:
	bsplits =  b.strip().split()
	# bimfile format: chr,=rsID, genPos, realPos
	chrom = int(bsplits[0])
	if (chrom in all_pos):
		all_pos[chrom].append((bsplits[1], int(bsplits[3])))
	else:
		all_pos[chrom] = [ (bsplits[1], int(bsplits[3]))]


## Then, build a dictionary of lists of tuples of chromosomal ranges for each individual to mask
## Dictionary indexed by chromosome 
print "Building masking dictionary..."
all_inds=list()
ind_keep = dict()
indfile = open(args.inds)
for i in indfile:
	isp = i.strip().split()
	all_inds.append( (isp[0],isp[1]))
	# Read in the bedfile
	if (args.diffNames):
		thisbed = open(isp[2], "r")
	else:
		thisbed = open(isp[1] + ".bed", "r")
	
	ind_keep[isp[1]] = dict()
	for p in thisbed:
		bedline = p.strip().split()
		chrom=int(bedline[0])
		if (bedline[3] in keep):
			if (chrom in ind_keep[isp[1]]):
				ind_keep[isp[1]][chrom].append( (int(bedline[1]), int(bedline[2])))
			else:
				ind_keep[isp[1]][chrom] = [ (int(bedline[1]), int(bedline[2])) ]


## Then, read through the bimfile line by line
## Make a list of rsids, position, and chromosome

## For each individual, for each rsid, check if it is in the correct range
## If it is, then keep it, if not, write to cluster file:
## for ind1:
## rsidxxx C1
## rsidxxx C1
print "Writing cluster files..."
## open a cluster file
myclust = 1
zerocluster = open("zerocluster.txt", "w")
## Read through each individual
for myind in all_inds:
	thisind = myind[1]
	clusterfile = open(thisind + "cluster.txt", "w")
	
	# Let's go through each chromosome
	for thischrom in all_pos:
		for record in all_pos[thischrom]:
			keep = False
			# Now we need to determine if this region is to be kept
			if thischrom in ind_keep[thisind]:
				for myinterval in ind_keep[thisind][thischrom]:
					if myinterval[0] < record[1] < myinterval[1]:
						keep = True
						break
			if not keep:
				clusterfile.write(record[0] +"\t" + "C" + str(myclust)+"\n")
	zerocluster.write(myind[0] + "\t" + thisind + "\tC" + str(myclust) +"\n")
	myclust += 1
	clusterfile.close()
zerocluster.close()
## for ind 2:
## rsidxxx C2
## rsidxxx C2
## and so on

## Need to generate the zero cluster file as well:
# ind1 ind1 C1
# ind2 ind2 C2

## Now that we have generated all the cluster files and individuals, run plink
print "Running plink missing..."
plinkhead = re.sub(r".bim$", "", args.bim)
thisplink =plinkhead
clust = 1
for myind in all_inds:
	thisout = plinkhead + "_C" + str(clust)
	clustfile = myind[1] + "cluster.txt"
	os.system("plink --bfile " + thisplink + " --zero-cluster " + clustfile + " --within zerocluster.txt --make-bed --out " + thisout)
	thisplink=thisout
	clust+=1



