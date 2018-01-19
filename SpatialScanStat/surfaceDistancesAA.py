### Program to compute surface distances between amino acids
### By J. Homburger
### September 2, 2015

# Import lots of packages
import argparse
import networkx as nx
import numpy as np
import sys
from scipy.spatial import cKDTree
from scipy.spatial.distance import euclidean

### Parse arguments ###
# Need to know:
# 1. PDB File that has been converted into a tab or space delimited table, this means all fields are separated by spaces
# 2. Vertex file output from MSMS
# 3. Faces file output from MSMS
# 4. Output file name

parser = argparse.ArgumentParser()

parser.add_argument("-p", help="pdb file name", required=True)
parser.add_argument("-v", help="MSMS vertex file", required=True)
parser.add_argument("-e", help="MSMS faces file", required=True)
parser.add_argument("-o", help="Output file name", required=True)
parser.add_argument("-d", help="Amino acid depth file name")
parser.add_argument("-start", help="Amino Acid Residue to Start at", default=1, type=int)
parser.add_argument("-finish", help="Amino Acid Residue to finish at", default=10000, type=int)
parser.add_argument("-center_search", help="Radius of search for best point to represent each amino acid", default=3, type=int)
parser.add_argument("-center_nearby", help="Number of points to search to represent each amino acid", default=100, type=int)
parser.add_argument("-chain", default="A", type=str)

args = parser.parse_args()

### Arguments in detail ###

## -p: pdb file name (needs to be space or tab delimited, many pdb files are not originally in this format
## -v: Vertex file output from the MSMS program
## -e: Faces file output from the MSMS program
## -o: name of the output file containing surface distances between the amino acids
## -d: name of the output file containing depths of each amino acid
## -start: amino acid residue to begin on, defaults to 1
## -finish: amino acid residue to finish with, defaults to 10000 or the highest amino acid
## Note: start/finish parameters can be used to easily parallelize the analysis (recommended)

## -center_search: when determining the closest surface point to the amino acid, this is the radius that potential points can be in
## -center_nearby: when determining the closest surface points to the amino acid, this is the maximum number of points to check
## The two commands above change the search heuristics. By making these assumptions we 
## speed up runtime significantly with little loss in accuracy.

## -chain: default amino acid chain to search. Defaults to chain A
## Note: can only do 1 chain at a time due to residue number conflicts

####### Main Script #######

# Read in the vertices

vfile = open(args.v, "r")

# Read past the 3 header lines
vfile.readline()
vfile.readline()
vparams = vfile.readline().strip().split()

vertices = np.empty([int(vparams[0]),3])

i = 0
print "Building the array of vertices..."
for vline in vfile:
	vsplit = vline.strip().split()
	## Positions [0, 1, 2] are the x,y,z coordinates
	vertices[i] = np.array( map(float,(vsplit[0], vsplit[1], vsplit[2])))
	i += 1
vfile.close()
### Read in the edges


efile = open(args.e, "r")

# Represent the edges as a set of tuples
# Remember, edges are 1-indexed, so need to knock 1 off the 
# values when reading them in to conform to the vertex array
edges = set()

# Read the first three lines, these are header rows
efile.readline()
efile.readline()
efile.readline()

print "Identifying the set of edges..."
for eline in efile:
	## Parse the line
	esplits = eline.strip().split()
	# Sort to ensure that the lowest is always first
	tri = map(int, esplits[0:3])
	tri.sort()
	# Remember to knock off 1 so that the indices match
	edges.add((tri[0] - 1, tri[1] - 1))
	edges.add((tri[0] - 1, tri[2] - 1))
	edges.add((tri[1] - 1, tri[2] - 1))

efile.close()
### Read in the pdb file ###

pdbfile = open(args.p, "r")

# List of tuples with residue number as key and list of np array coordinates for each atom in a residue
allatoms = dict()

# Now this I don't think has a standard number of lines on the top
print "Reading in the PDB file..."
for pline in pdbfile:
	# Check to make sure we only read ATOM lines
	if not pline.startswith('ATOM'):
		continue
		
	## So, the order is as follows
	## Atom, AtomNum, AtomType, ResidueType, Chain, ResidueNum, X, Y, Z, Unk, Beta, Unk
	psplits = pline.strip().split()
	
	## Need to ensure its the correct chain
	if not psplits[4] == args.chain:
		continue
	## Get xyz coordinates
	xyz = np.array(map(float, [psplits[6], psplits[7], psplits[8]]))
	if int(psplits[5]) in allatoms:
		allatoms[int(psplits[5])].append( xyz)
	else:
		allatoms[int(psplits[5])] = [ xyz ]



pdbfile.close()

### Now assign the atoms to the closest vertex

# Build the tree structure 
disttree = cKDTree(vertices)

# Dictionary mapping a residue to the closest vertex
closest_point = dict()
# Dictionary holding a residues minimum distance to the protein surface 
aa_depth = dict()

print "Calculating distance to surface for each amino acid residue..."
# This also assigns each amino acid residue to the "closest" surface point
for z in allatoms.keys():
	k = allatoms[z]
	# Does it already exist in the data frame?
	
	# calculate the point that has the minimum sum of distances to all of the non-hydrogen atoms in the molecule
	
	minx = None
	minvalue = np.inf
	checkpoints = set()
	for thispoint in k:
		thislookup = disttree.query_ball_point(thispoint, r=args.center_search)
		checkpoints.update(thislookup)

	if len(checkpoints) == 0:
		lookups = disttree.query(k[0], k=args.center_nearby)
		checkpoints = lookups[1]
	
	#print len(checkpoints)
	#print z
	thisdepth = np.inf
	for x in checkpoints:
		totaldist = 0.0
		
		for thispoint in k:
			totaldist += euclidean(thispoint, vertices[x])
		if totaldist < minvalue:
			minx= x
			minvalue = totaldist
			thisdepth = totaldist/len(k)
	aa_depth[z] = thisdepth
	closest_point[z] = minx
	#print z


#print closest_point

### Build the graph

surfaceG = nx.Graph()
print "Building the weighted graph structure..."
# Add the nodes
surfaceG.add_nodes_from(range(0, vertices.shape[0]))

print "Adding edges to the graph..."
# Add the edges
for e in edges:
	# For each edge, get the two vertex indices, and calculate the distance
	thisdist = euclidean(vertices[e[0]], vertices[e[1]])
	surfaceG.add_edge(e[0], e[1], weight=thisdist)

print "Graph generated, calculating distances..."
# We need to do a "bottom triangle" distance calculation

# Get a list of all the residue numbers

residues = closest_point.keys()
residues.sort()

# Now let's go through the entire list
print "Calculating distance between all pairs of amino acids..."

## Need a heuristic function for the a_star search, going to use the Euclidean distance between the two nodes

def astar_heur_euc(u,v):
	return(euclidean(vertices[u], vertices[v]))

surface_distances = dict()
start_aa = max(args.start, 1)
end_aa = min(args.finish+1, max(residues)+1)
for i in range(start_aa,end_aa):
#for i in range(1,2):
	sys.stdout.write("Processing residue number: %s   \r" % (i) )
	sys.stdout.flush()
	if not i in closest_point:
		continue
	for j in range(i, max(residues)+1):
#	for j in range(i, 100):
		if not j in closest_point:
			continue
	#	sys.stdout.write("Processing residue number: %s  compared to %s \r" % (i,j) )
	#	sys.stdout.flush()
		try:
			thislength = nx.astar_path_length(surfaceG, closest_point[i], closest_point[j], weight='weight', heuristic=astar_heur_euc)
		except nx.NetworkXNoPath:
			thislength = float("inf")
		surface_distances[ (i,j)] = thislength
sys.stdout.write("\n\n")
print "Writing distance output to file"

outfile = open(args.o, "w")

# Output surface distances to file
for i in range(start_aa,end_aa):
#for i in range(1,2):
	for j in range(i, max(residues)+1):
#	for j in range(i, 100):
		if not (i,j) in surface_distances:
			continue
		outfile.write( "\t".join(map(str, [i, j, surface_distances[ (i,j)]])) + "\n" )	
outfile.close()


# Output amino acid depths, only done if filename is specified in arguments
if (args.d):
	print "Writing calculated amino acid depths to file"
	surfout = open(args.d, "w")
	for i in residues:
		if not i in aa_depth:
			continue
		surfout.write( str(i) + "\t" + str(aa_depth[i]) + "\n")





















