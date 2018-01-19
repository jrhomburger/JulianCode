### Second IBD Pipeline 

## Goal: Take phased data, run GERMLINE, and output a list of all the IBD matches

## Inputs:
# 1. 22 Phased SHAPEIT Style files (potentially gzipped)
# 2. BIM file from analysis
# 3. FAM file from analysis
# 4. SNP Locations/Genetic map information
# This version should only be run one chromosome at a time

import re
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-s', '--shpfile', help='SHAPEIT output file name')
parser.add_argument('-f', '--famfile', help='plink FAM file')
parser.add_argument('-m', '--bimfile', help='plink BIM file')
parser.add_argument('-l', '--locfile', help='SNP reference map file')
parser.add_argument('-o', '--outfile', help="GERMLINE outfile prefix", default="gline")
parser.add_argument('-maf', '--maf', default=0.0, type=float, help='MAF cutoff, SNPs below this MAF will be removed')
parser.add_argument('-minl', '--minl', default=5, type=float, help='Minimum length of IBD matches to be reported')
parser.add_argument('--gpath', default="germline", type=str, help="Path or call for Germline function")
parser.add_argument('-haploid')


commonvariants = set()

def makeSNPMap(bimname, referencemap):
	"""
	Function that given prefixes of the file of locus information and a reference genetic map.
	Creates a dictionary indexed by position of snp genetic positions
	Returns this dictionary
	Please only do one chromosome at a time...
	"""
	
	snpmap = dict()

		# Open chromosome map files
		
	bimfile = open(bimname, "r")
	
	# Open map files, swap 22 with current chromosome
	
	mapfile = open(referencemap, "r")
	
	# Initialize variables 
	previousCM = 0
	previousPos = 0
	i=0
	bimline = bimfile.readline().strip().split() # Pos 1 is rsID, Pos 3 is location
	for mapline in mapfile:
		if len(bimline) == 0:
			break		
		if i==0:
			i+=1
			continue
		mapline = mapline.strip().split()
		# Three cases: 1. SNP pos gt map pos
		while int(bimline[3]) < int(mapline[0]): # This means that the BIM file is behind the map file, so need to write output here with the interopolation
		# of the previous values
			diffCM = float(mapline[2]) - float(previousCM)
			diffpos = float(mapline[0]) - float(previousPos)
			multi = (float(bimline[3]) - float(previousPos))/diffpos
			cmout = multi*diffCM + float(previousCM)
			if cmout < 0: # this should not happen so if it does dump data and quit
				print i
				print cmout
				print diffCM
				print diffpos
				print previousCM
				print previousPos
				print bimline
				print mapline
				exit()

			snpmap[bimline[3]] = cmout
			bimline = bimfile.readline().strip().split()
			if len(bimline) == 0:
				break		
		if len(bimline) ==0:
			break
		if bimline[3] == mapline[0]: # write out genetic position
			snpmap[bimline[3]] = cmout
			bimline = bimfile.readline().strip().split()

		#if bimline[3] > mapline[0]: # read next line in the map file
		#	previousCM = mapline[2]
		#	previousPos = mapline[0]
		#	continue
		# Hits this and continues if bimline is above mapline
		previousCM = mapline[2]
		previousPos = mapline[0]
		i += 1
	bimfile.close()
	mapfile.close()
	return(snpmap)

def maf(alleles):
	tester = float(alleles.count(alleles[0]))/len(alleles)
	return(min([tester, 1.0-tester]))
	
def germline_ibd(pedfile, mapfile, outfile, haploid=1, bits=128, min_m = 3, err_hom = 2, err_het = 1, call="germline "):
	"""
	parameters:
	in_prefix: prefix for the .ped and .map input files
	haploid: flag to run haploid version of GERMLINE
	bits: bits used for seeding the match
	min_m: minimum cM length of a match
	err_hom: homozygous error tolerance in a match
	err_het: heterozygoues error tolerance in a match
	gline_path: path to the gline.sh file
	
	This function assumes the .ped files is phased and the .map file contains genetic positions in cM
	Also going to assume that the germline program is in the current path
	I also recommend against including ; rm -rf *; semicolon in the file input name	
	
	Germline will output the file inprefix_gout.match which is a table of matches and inprefix_gout.log 
	
	Note: GERMLINE works best when only a single chromosome is run at a time
	"""
	
	# Write filenames to the run file
	runfile = open("out.run", "w")
	runfile.write("1\n")
	runfile.write(mapfile + "\n")
	runfile.write(pedfile + "\n")
	runfile.write(outfile + "\n")
	runfile.close()
	
	# Call program and add files
	
	# Update the parameters
	if haploid:
		call = call + " -min_m " + str(min_m) + " -haploid -bits " + str(bits) + " -err_hom " + str(err_hom) + " -err_het " + str(err_het) + " < out.run"
	else:
		call = call + " -min_m " + str(min_m) + " -bits " + str(bits) + " -err_hom " + str(err_hom) + " -err_het " + str(err_het) + " < out.run"
	print(call)
	os.system(call)
	os.remove("out.run")


### Part 1: Convert SHAPEIT Files into compatible Plink .ped files by chromosome

args = parser.parse_args()

shpfile = args.shpfile
bimfile = args.bimfile
famfile = args.famfile
locfile = args.locfile
minl = args.minl

# First make SNP location dictionary - to add 

snpmap = makeSNPMap(bimfile, locfile)

# Next, read in the alleles from the HAPs file

print 'reading annotation information from',bimfile
print 'line:'
mapdict = {}
last_pos = 0

for i, line in enumerate(file(bimfile)):
	print '\b\b\b\b\b\b\b\b\b\b%.8g' % (i),
	sys.stdout.flush()
	line = line.strip().split()
	mapdict[line[1]] = line[:4]
	#if int(line[0]) != cur_chr:
	#	chrdict[cur_chr] = (last_pos, i -1)
	#	last_pos = i
	#	cur_chr = int(line[0])

print '\nreading pedigree/phenotype information from',famfile
famdata = [line.strip().split() for line in file(famfile).readlines()]

alleles = []
snplist = []
indlist = []
keep = []
last_chr = 0
cur_start = 0
below_maf = 0
j=0
print 'reading in genotypes...'
for i,line in enumerate(file(shpfile)):
	#print '\b\b\b\b\b\b\b\b%.6g' % (i),
	#sys.stdout.flush()
	line = line.strip().split()
	thismaf = maf(line[5:])
	if thismaf < args.maf:
		below_maf += 1
		continue
	j+=1
	snplist.append(line[1])
	allele_ones="".join(line[5:])
	allele_fix1 = re.sub(r"0", line[3], allele_ones)
	allele_fix2 = re.sub(r"1", line[4], allele_fix2)
	alleles.append(list(allele_fix2))
	
# Output only variants above a certain MAF


print "\nI removed " + str(below_maf) + " SNPs due to the MAF cutoff"
print "\nadding SNP map positions to the SNP information..."

for s in snplist:
	sinfo = mapdict[s]
	if int(sinfo[3]) in snp_cm:
		mapdict[s][2] = snpmap[int(sinfo[3]) ]
	else:
		print "WARNING: I couldn't find a genetic position for the following SNP:\n"
		print (int(sinfo[0]), int(sinfo[3]))

print '\ntransposing alleles to individuals-in-rows format...'
alleles = zip(*alleles)

outfilename=args.outfile+".ped"
outmapname= args.outfile+".map"

outfile = file(outfilename,'w')
outmap = file(outmapname,'w')


print '\n\nwriting output to',outfilename
# assert twice as many alleles as fam data - will only work reliably for autosomes
# But good idea to check this just to make sure...
assert len(alleles) == 2*len(famdata)
assert 2*len(famdata) == len(indlist)
nsnps = len(snplist)
for i in range(len(famdata)):
   print '\b\b\b\b\b\b\b\b%.5g' % (i),
   thisname = re.sub(r'_[AB]', '', indlist[i])
   sys.stdout.flush()
   newline = list()
   for k in range(0, len(alleles)):
	   #print k
	   a1 = alleles[2*i][k]
	   a2 = alleles[2*i + 1][k]
	   newline.append(a1)
	   newline.append(a2)
   outfile.write('%s\t%s\n' % ('\t'.join(famdata[i]),'\t'.join(newline)))

print '\n\nwriting map file output to',outmapname
for i in range(len(snplist)):
   #print '\b\b\b\b\b\b\b\b%.6g' % (i),
   #sys.stdout.flush()
   if getchr[snplist[i]] == c:
	   outmap.write('%s\n' % ('\t'.join(mapdict[snplist[i]])))
# Close your goddamn filehandles
outfile.close()
outmap.close()

### Now run Germline on this individual chromosome

outname = args.outfile+".gline"
hap=0
if (args.haploid):
	hap=1

germline_ibd(outfilename, outmapname, outname, min_m = str(args.minl), haploid = hap, call=args.gpath)


