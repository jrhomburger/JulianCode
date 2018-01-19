### Integrated script for running IBD analyses 

## Goal: Take phased data, run GERMLINE, and output a list of all the IBD matches

## Inputs:
# 1. Phased beagle file
# 2. BIM file from analysis
# 3. FAM file from analysis
# 4. SNP Locations/Genetic map information

## Outputs:
# GERMLINE file of IBD matches for all chromosomes

## Thanks to C. Gignoux for some of the original code in here. 

import sys
import argparse
import re
import subprocess
import math
parser = argparse.ArgumentParser()

parser.add_argument('-b', '--bglfile', help='BEAGLE output filename')
parser.add_argument('-f', '--famfile', help='plink FAM file')
parser.add_argument('-m', '--mapfile', help='plink MAP file')
parser.add_argument('-l', '--locfile', help='SNP locations file prefix, this will be XXXX_chrY.snp_locations')
parser.add_argument('-refmap', '--refmap', help="Reference mapfile, used to create .snp_locations file if those aren't present")
parser.add_argument('-maf', '--maf', default=0.0, type=float, help='MAF cutoff, SNPs below this MAF will be removed')
parser.add_argument('-minl', '--minl', default=5, type=float, help='Minimum length of IBD matches to be reported')
parser.add_argument('-haploid')

commonvariants = set()

def maf(alleles):
	tester = float(alleles.count(alleles[0]))/len(alleles)
	return(min([tester, 1.0-tester]))
	
def germline_ibd(pedfile, mapfile, outfile, haploid=1, bits=128, min_m = 3, err_hom = 2, err_het = 1):
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
	"""
	
	# Write filenames to the run file
	runfile = open("out.run", "w")
	runfile.write("1\n")
	runfile.write(mapfile + "\n")
	runfile.write(pedfile + "\n")
	runfile.write(outfile + "\n")
	runfile.close()
	
	# Call program and add files
	call = "germline "
	
	# Update the parameters
	if haploid:
		call = call + " -min_m " + str(min_m) + " -haploid -bits " + str(bits) + " -err_hom " + str(err_hom) + " -err_het " + str(err_het) + " < out.run"
	else:
		call = call + " -min_m " + str(min_m) + " -bits " + str(bits) + " -err_hom " + str(err_hom) + " -err_het " + str(err_het) + " < out.run"
	print(call)
	os.system(call)
	os.remove("out.run")

def makeSNPMap(snpfile, referencemap):
	"""
	Function that given a file of locus information and a reference genetic map.
	Creates a cM delimited mapfile for the given loci.
	Returns the name of the mapfile.
	"""
	bimfile = open(snpfile, "r") # open the input file
	mapfile = open(referencemap, "r")
	outfilename = re.sub(r'\.bim', '.markerpos', snpfile)
	posfilename = re.sub(r'\.bim', '.snp_locations', snpfile)
	outfile = open(outfilename, "w")
	posfile = open(posfilename, "w")
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

			outfile.write( str(cmout) +"\n")
			posfile.write( str(bimline[3]) + "\t" + str(cmout) + "\n")
			bimline = bimfile.readline().strip().split()
			if len(bimline) == 0:
				break		
		if len(bimline) ==0:
			break
		if bimline[3] == mapline[0]: # write out genetic position
			outfile.write( mapline[2]+ "\n")
			posfile.write( str(bimline[3]) + "\t" + mapline[2] + "\n")
			bimline = bimfile.readline().strip().split()
	
		#if bimline[3] > mapline[0]: # read next line in the map file
		#	previousCM = mapline[2]
		#	previousPos = mapline[0]
		#	continue
		# Hits this and continues if bimline is above mapline
		previousCM = mapline[2]
		previousPos = mapline[0]
		i += 1
	outfile.close()
	return(outfile.name)


### Part 1: Convert BEAGLE File into compatible Plink Files by chromosome

args = parser.parse_args()

bglfile = args.bglfile
mapfile = args.mapfile
famfile = args.famfile
minl = args.minl

print 'reading from bgl file',bglfile
#print 'will write output to %s and %s' % (outfilename,outmapname)

print 'reading annotation information from',mapfile
print 'line:'
mapdict = {}
chrdict = {}
getchr = {}
cur_chr = 1
last_pos = 0

for i, line in enumerate(file(mapfile)):
	print '\b\b\b\b\b\b\b\b\b\b%.8g' % (i),
	sys.stdout.flush()
	line = line.strip().split()
	mapdict[line[1]] = line[:4]
	getchr[line[1]] = int(line[0])
	#if int(line[0]) != cur_chr:
	#	chrdict[cur_chr] = (last_pos, i -1)
	#	last_pos = i
	#	cur_chr = int(line[0])

chrdict[cur_chr] = (last_pos, i)
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
for i,line in enumerate(file(bglfile)):
	#print '\b\b\b\b\b\b\b\b%.6g' % (i),
	#sys.stdout.flush()
	
	#only care about lines that have markers - fam data comes from famfile
	if line.startswith('I'):
		indlist.extend(line.strip().split()[2:])
	
	
	if line.startswith('M'):
			line = line.strip().split()
			thismaf = maf(line[2:])
			if thismaf < args.maf:
				below_maf += 1
				continue
			j+=1
			snplist.append(line[1])
			alleles.append(line[2:])

			#3 Need to build chr dictionary based on this, so i starts at +1 
			if getchr[line[1]] != last_chr:
				chrdict[last_chr] = (cur_start-1, j - 2) 
				#print last_chr
				#print getchr[line[1]]
				#print i
				#print cur_start
				last_chr = getchr[line[1]]
				#print last_chr
				#sys.exit()
				cur_start = j
				
## Dictionary holds the list boundaries of each chromosome 		
chrdict[last_chr] = (cur_start-1, j - 1)			
print chrdict
print "\nI removed " + str(below_maf) + " SNPs due to the MAF cutoff"
print "\nadding SNP map positions to the SNP information..."
## Add positions to the SNP list

if (!args.locfile):


	locfile = 
else:
	locfile = args.locfile


snp_cm = dict() # stores the genetic positions with a tuple of (chromosome, snp)
for c in range(1,23):
	thisloc = open(locfile + "_chr" + str(c) + ".snp_locations")
	for f in thisloc:
		fsplit = f.strip().split()
		snp_cm[ (c, int(fsplit[0]))] = fsplit[1]
		
## This adds the positions from the dictionary into the current SNP list
## Added genetic positions will be printed in output
for s in snplist:
	sinfo = mapdict[s]
	if (int(sinfo[0]), int(sinfo[3])) in snp_cm:
		mapdict[s][2] = snp_cm[ (int(sinfo[0]), int(sinfo[3])) ]
	else:
		print "WARNING: I couldn't find a genetic position for the following SNP:\n"
		print (int(sinfo[0]), int(sinfo[3]))

print '\ntransposing alleles to individuals-in-rows format...'
alleles = zip(*alleles)

allpeds = []
allmaps = []

# For each chromosome, make a new output file
for c in range(1,23):

	outfilename = bglfile.replace('.beagle','_chr' + str(c) + '.ped')
	outmapname = bglfile.replace('.beagle','_chr' + str(c) + '.map')
	allpeds.append(outfilename)
	allmaps.append(outmapname)
	
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
		chrstart = chrdict[c][0]
		chrend = chrdict[c][1]
		for k in range(chrstart, chrend+1):
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

### Part 2: Run GERMLINE and concatenate the output


all_outfiles = []



print allpeds
print allmaps
#sys.exit()
# For each ped/map combo, run Germline
for z in range(len(allpeds)):
	thisped=allpeds[z]
	thismap=allmaps[z]
	outname = thisped.replace(".ped", "_diploid_gline")
	hap=0
	if (args.haploid):
		hap=1
	
	germline_ibd(thisped, thismap, outname, min_m = str(args.minl), haploid = hap)
	all_outfiles.append(outname + ".match")



## Concatenate all of the output files together
gline_outfilename = bglfile.replace('.beagle','.gline')
cat_call = "cat " + " ".join(all_outfiles) + " > " + gline_outfilename
subprocess.call(cat_call, shell=True)





	
	
	
	



