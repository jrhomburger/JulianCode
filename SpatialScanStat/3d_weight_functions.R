### #d weighting functions
### These functions are used to generate colors that are smoothed
### across a protein model. So one can view some distribution of 
### data across the spatial gradient of a protein in a program
### such as Pymol

## Nadaraya-Watson Kernel smoothing average: 
## Ie. Local average with values weighted by distance to point
## The trick is to get the right distance kernel

## Epanichnikov kernel:
## 3/4 * ( 1 - u^2) * indicator if u is less than lambda 

# let u be the euclidean norm divided by lambda scale parameter

euclid <- function(x1, x2) {
  # Calculates euclidean distances between points
  return(sqrt(sum((x1-x2)^2)))
}

eKernel <- function(dist, lambda) {
  # Should return the correct epanichnikov weights
  u = dist/lambda
  return( (3/4)*(1 - u^2)*(abs(u) < 1) )
}

genPosList <- function(residues, pdb_index, pdb_coords) {
  # given the pdb frame, returns a list indexed by residue that returns the 3d point positions
  # PDB Index and PDB Coords need to be in the same order - ie. element 1 of pdb index corresponds
  # to row 1 of the PDB Coords.
  # PDB Coords must also be a matrix, even if it is only one dimension
  indices <- residues
  positions = list()
  for (i in indices) {
    if (!i%in%pdb_index) {
      positions[[paste("pos",as.character(i), sep="")]] <- c(NA,NA,NA)
    } else {
      j = which(pdb_index==i)[[1]]
      positions[[paste("pos",as.character(i), sep="")]] <- as.vector(pdb_coords[j,])
    }
  }
  return(positions)
}

# Let's functionalize the kernel smoothing method for future ease

KernelEstimate3D <- function(phenotype_values, phenotype_indexes, all_positions, estimate_indexes, lambda=10) {
  ### Input:
  
  # phenotype_values: vector of observed phenotype values
  # phenotype indexes: list of the 'index' of each of the phenotypes, for our protein example this is the amino acid residue
  # all_positions: list, indexed by the indexes, of the 3d (or however many dimensions) positions of all the residues. NA represents missing position info
  # estimate_indexes: list of the indexes to get kernel weighted averages at
  # lambda: tuning parameter for epanovichnikov kernel 
  
  output_vector <- rep(NA, length(estimate_indexes)) 
  weight_vector = rep(NA, length(estimate_indexes))
  for (i in 1:length(estimate_indexes)) {
    residue = estimate_indexes[i]
    pos = all_positions[[paste("pos",residue,sep="")]]
    if (is.na(pos[1])) {
      output_vector[i] = NA
    }
    curtop = 0
    curbot = 0
    for (j in 1:length(phenotype_values)) {
      ## Calculate the euclid distance between the two
      if (is.na(phenotype_values[j])) {
        next
      }
      if (is.na(phenotype_indexes[j]) | !paste("pos",phenotype_indexes[j], sep="")%in%names(all_positions)) {
        next
      } # Skip if its out of range:jj <- 301
      thisdist <- euclid(pos, all_positions[[paste("pos",phenotype_indexes[j], sep="")]])
      if (is.na(thisdist) | thisdist/lambda > 1) { # Out of the kernel range
        next
      }
      thiskernel = eKernel(thisdist, lambda)
      curtop = curtop + thiskernel*phenotype_values[j]
      curbot = curbot + thiskernel
      
    }
    
    output_vector[i] = curtop/curbot
    weight_vector[i] = curbot
    if (i %%10 == 0) {
      print(i)
    }
    
  }
  
  return(data.frame(index=estimate_indexes, value=output_vector, weight=weight_vector))
  
}

kernelColors <- function(weighted_values, low_cutoff = 0, high_cutoff = 9999, num_levels = 20, high_col = "red", low_col = "blue", exclude_above=10000, exclude_below=0, exclude_indices=c(NA)) {
  # weighted values is a list of the weighted values. the other parameters are:
  # low_cutoff, high_cutoff = cutoffs for out of range values in the weighted list. Use if having regions of high variance
  # num_levels = number of color levels to cut the data into, don't know what the upper limit is but I generally use 20
  # high_col, low_col = names of colors used for high values and low values respectively
  # exclude above and exclude below: list indexes to exclude above and below, useful for dealing with edge effects
  # exclude indices: indices to exclude
  age_colors_fxn <- colorRampPalette(c(low_col, high_col))
  wt2 = weighted_values
  wt2[wt2 > high_cutoff | wt2 < low_cutoff] <- NA
  wt2[is.nan(wt2)] <- NA
  if (exclude_below>0) {
    wt2[1:exclude_below] <- NA
  }
  if (exclude_above <= length(wt2)) {
    wt2[exclude_above:length(wt2)] <- NA
  }
  
  age_colors <- age_colors_fxn(num_levels)[as.numeric(cut(wt2, breaks=num_levels))]
  
  age_colors[is.na(age_colors)] <- "#A8A8A8"
  
  return(age_colors)
}

## Read in a PDB file from the protein data bank
pdb_file <- read.table("")
pdb_ca <- subset(pdb_file, V3=="CA")
names(pdb_ca) <- c("type", "atomID", "atomType", "AA", "Chain", "ResidueNo", "X", "Y", "Z", "Unk", "Temp", "Unk2")

some_data_file <- read.table("") # Some data file with values of interest and
# protein positions, expects names value and amino acid number  

protein_positions <- genPosList(1:841, pdb_file$ResidueNo, matrix( c(pdb_file$X,  pdb_file$Y,  pdb_file$Z), ncol=3))
smoothing_param = 30
KernelEstimate3D(some_data_file$value, some_data_file$amino_acid_number, protein_positions, 1:max(pdb_ca$ResidueNo), lambda=smoothing_param)

