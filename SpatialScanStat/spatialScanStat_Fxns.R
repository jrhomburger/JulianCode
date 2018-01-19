### These are the functions used to calculate the spatial scan statistic over protein data
### The spatial scan statistic is developed in Kuldorff et. al. 1997. Here we apply it to
### disease true/false states for each genomic mutation. We apply this across the protein's
### three dimensional structure.

### This also contains code for performing this scan across the surface of a protein and
### for munging some data tables into the correct format.


lik_sss_pseudo = function(yw, yg, nw, ng  ) {
  ## Function calculates the likelihood given some counts, define:
  ## yw = number of affected rare variants in region
  ## yg = number of affected rare variants in whole gene/region
  ## nw = number of rare variants in the window
  ## ng = number of rare variants in whole gene/region
  
  ## Determine calculables: 
  
  # Is the pseudocounting ok?
  pw_hat = (yw+1)/(nw+2) 
  qw_hat = (yg - yw + 1)/(ng - nw + 2)
  rg_hat = (yg+1)/(ng+2)
  
  ## Calculate the log of the likelihood ratio statistic for the region
  if (pw_hat > qw_hat) {
    llrw = yw * log(pw_hat/rg_hat) + (nw - yw)*log ( (1-pw_hat)/(1-rg_hat) ) + (yg - yw ) * log(qw_hat/rg_hat) + (ng - nw - (yg-yw) )*log( (1-qw_hat)/(1-rg_hat))
  } else {
    llrw = 0
  }
  
  return(llrw)  
}

lik_sss = function(yw, yg, nw, ng  ) {
  ## Function calculates the likelihood given some counts, define:
  ## yw = number of affected rare variants in region
  ## yg = number of affected rare variants in whole gene/region
  ## nw = number of rare variants in the window
  ## ng = number of rare variants in whole gene/region
  
  ## Determine calculables: 
  
  # Is the pseudocounting ok?
  if (nw!=0) {
    pw_hat = (yw)/(nw) 
  } else {
    return(0)
  }
  if ((ng-nw)==0) {
    ## All sites outside the window are affected
    return(0)
  }
  
  
  qw_hat = (yg - yw)/(ng - nw)
  rg_hat = (yg)/(ng)
  
  ## Calculate the log of the likelihood ratio statistic for the region
  if (pw_hat > qw_hat & (1-pw_hat) == 0) {
    llrw = yw * log(pw_hat/rg_hat) + 0 + (yg - yw ) * log(qw_hat/rg_hat) + (ng - nw - (yg-yw) )*log( (1-qw_hat)/(1-rg_hat))
  } else if (pw_hat > qw_hat) {
    #if( pw_hat > qw_hat) {
    llrw = yw * log(pw_hat/rg_hat) + (nw - yw)*log ( (1-pw_hat)/(1-rg_hat) ) + (yg - yw ) * log(qw_hat/rg_hat) + (ng - nw - (yg-yw) )*log( (1-qw_hat)/(1-rg_hat))
  } else {
    llrw = 0
  }
  
  return(llrw)  
}

eucdist = function(a,b) {
  ## Returns the Euclidean distance between two points
  ## Doesn't sanitize for dimensions so be careful
  return ( sqrt( sum((a-b)^2)))
}

## Window Functions
eucWinLik = function( affected, position, window_size, window_center) {
  ### Take in data with at least the following inputs:
  ## 1. Affected status - assumes this is a Boolean True/False or 1/0 vector
  ## 2. Position (in whatever dimension - needs to be coerced to a matrix) - also same
  ## order as the affected status
  ## 3. Window center - position (in same dimensions as the position) that is the center
  ## of the window to look in
  ## 4. Distance - maximum Euclidean distance to look in the window from
  
  
  ## First, generate subset of within the window
  inWindow = unname(apply(position, 1, eucdist, b=window_center) <= window_size)
  ## Let's get the count of the mutations inside the window
  yw = sum(inWindow & affected)
  nw = sum(inWindow)
  
  ## Total affected and mutations outside of window
  yg = sum(affected)
  ng = length(affected)
  
  thislik = lik_sss(yw, yg, nw,ng)
  
  # Returns likelihood and size of window (in AAs)
  return(list(lik=thislik, windowSize=sum(inWindow)))
}

eucWinLik_sparse = function( affected, position, window_size, window_center, pseudo_count=F) {
  ### Take in data with at least the following inputs:
  ## 1. Affected status - assumes this is a Boolean True/False or 1/0 vector
  ## 2. Position (in whatever dimension - needs to be coerced to a matrix) - also same
  ## order as the affected status
  ## 3. Window center - position (in same dimensions as the position) that is the center
  ## of the window to look in
  ## 4. Distance - maximum Euclidean distance to look in the window from
  
  
  ## First, generate subset of within the window
  inWindow = unname(apply(position, 1, eucdist, b=window_center) <= window_size)
  ## Let's get the count of the mutations inside the window
  yw = sum(inWindow & affected)
  nw = sum(inWindow)
  
  ## Total affected and mutations outside of window
  yg = sum(affected)
  ng = length(affected)
  if (pseudo_count) {
    thislik = lik_sss_pseudo(yw, yg, nw,ng)
  } else {
    thislik = lik_sss(yw, yg, nw,ng)
  }
  
  
  # Sparse functions return only the likelihood
  return(thislik)
}

## Single window_size functions
oneWinSize_sparse = function(affected, position, window_locations, window_size, window_func = eucWinLik_sparse, pseudo_count=F) {
  # Affected - T/F of affected status
  # Position - Matrix of positions
  # Window locations - List of the window locations
  # Window function - function to calcualte likelihood in windows, default is euclidean distance
  # Window size - size of window
  
  
  resout <- apply(window_locations, 1, window_func, affected=affected, position=position, window_size=window_size, pseudo_count=pseudo_count )
  
  ## Return the vector of the window test statistics
  return(resout)
}

oneWinSize = function(affected, position, window_locations, window_size, window_func = eucWinLik) {
  # Affected - T/F of affected status
  # Position - Matrix of positions
  # Window locations - List of the window locations
  # Window function - function to calcualte likelihood in windows, default is euclidean distance
  # Window size - size of window
  
  a = proc.time()
  resout = data.frame(t(sapply(window_locations, window_func, affected=affected, position=position, window_size=window_size )))
  
  ## Return the vector of the window test statistics
  return(list(test_stats=unlist(resout$lik), win_size = unlist(resout$windowSize)))
}

## Test functions
getMaxLik = function(affected, position, window_locations, window_sizes, window_func = eucWinLik) {
  ## Function calculates the maximum likelihood observed over the set of windows
  ## Returns a list of max stat, test stats, window locations, sizes and coverage
  ## See binomial model in Ionita-Laza et. al. 2012, AJHG and Kulldorff 1997 or any of Kulldorfs programs
  ## Inputs:
  # affected = T/F vector of affected status
  # position = vector of positions in the same order as the affected status
  # window_locations = locations to start a given window at
  # window_sizes = vecotr of sizes (in Angstroms) to test in spherical coordinates around each of the window locations
  # window_function = the specific window function to use, default uses Euclidean distances
  # pseudo_count = should we use pseudo-counts in the likelihood? default is False. 
  
  ### Note - the _sparse functions are better optimized and faster, they don't however return the window coverage
  
  allTestStats = NULL
  allWinSizes = NULL
  locs_for_func = as.list(data.frame(t(window_locations)))
  
  for (s in window_sizes) {
    
    thisiter <- oneWinSize(affected, position, locs_for_func, window_size = s, window_func = window_func)
    
    allTestStats = cbind(allTestStats, thisiter$test_stats)
    allWinSizes = cbind(allWinSizes, thisiter$win_size)
  }
  
  return( list(max_stat = max(allTestStats), test_stats = allTestStats, loc = window_locations, sizes = window_sizes, winCoverage = allWinSizes))
}

getMaxLik_sparse = function(affected, position, window_locations, window_sizes, window_func = eucWinLik_sparse, pseudo_count=F) {
  ## Function calculates the maximum likelihood observed over the set of windows
  ## Returns a list of max stat, test stats, window locations, sizes
  ## See binomial model in Ionita-Laza et. al. 2012, AJHG and Kulldorff 1997 or any of Kulldorfs programs
  ## Inputs:
  # affected = T/F vector of affected status
  # position = vector of positions in the same order as the affected status
  # window_locations = locations to start a given window at
  # window_sizes = vecotr of sizes (in Angstroms) to test in spherical coordinates around each of the window locations
  # window_function = the specific window function to use, default uses Euclidean distances
  # pseudo_count = should we use pseudo-counts in the likelihood? default is False. 
  
  allTestStats = matrix(NA, nrow=nrow(window_locations), ncol=length(window_sizes))
  
  for (s in 1:length(window_sizes)) {
    
    thislik <- oneWinSize_sparse(affected, position, window_locations, window_size = window_sizes[s], window_func = window_func,pseudo_count=F)
    
    allTestStats[,s] = thislik
  }
  
  return( list(max_stat = max(allTestStats), test_stats = allTestStats, loc = window_locations, sizes = window_sizes))
}

#3 Permutation functions
getPermLiks = function(affected, position, window_locations, window_sizes, nperms =100, window_func = eucWinLik) {
  ### Deprecated Function, use the getPermLiks_rep function instead
  ## Same inputs as the getMaxLik function, however this gets permtuation values
  #3 Can take a long time to run
  ## Better to use
  obs_liks = numeric(nperms)
  
  t = proc.time()
  for (j in 1:nperms) {
    this_lik = numeric(length(window_sizes))
    for (i in 1:length(window_sizes)) {
      this_lik[i] = max(oneWinSize(sample(affected, length(affected)), position, window_locations, window_size = window_sizes[i], window_func = window_func)$test_stats)
      
    }
    obs_liks[j] = max(this_lik)
    if (j%%10 == 0) {
      print(j)
      print(obs_liks[j])
    }
    
  }
  ltime = proc.time() - t
  print(ltime)
  return(obs_liks)
}

getPermLiksRep = function(affected, position, window_locations, window_sizes, nperms =100, window_func = eucWinLik_sparse) {
  ## Better optimized function for getting the likelihoods from the perumtations, still takes a bit of time though
  ## Inputs follow the inputs of getMaxLik
  r = proc.time()
  rliks = replicate(nperms, singlePermLiks(affected, position, window_locations, window_sizes, window_func =  window_func))
  rtime = proc.time() - r
  print(rtime)
  return(rliks) 
}

singlePermLiks = function(affected, position, window_locations, window_sizes, window_func = eucWinLik) {
  ## Helper function for the replicate function above
  ## Tried apply but it doesn't really spee dup the loop (only 4 items or so)
  this_lik = numeric(length(window_sizes))
  for (i in 1:length(window_sizes)) {
    this_lik[i] = max(oneWinSize_sparse(sample(affected, length(affected)), position, window_locations, window_size = window_sizes[i], window_func = window_func))
    
  }
  return(max(this_lik))
}

getAvgResidues = function(position, window_locations, window_sizes ) {
  ## Goal - for each window size, calculate the average number of residues that the window contains
  
  all_means = numeric(length(window_sizes))
  all_sds = numeric(length(window_sizes))
  for (s in 1:length(window_sizes)) {
    
    all_vals = numeric(length(window_locations))
    for ( j in 1:length(window_locations)) {
      all_vals[j] = sum(unname(apply(position, 1, eucdist, b=window_locations[j]) <= window_sizes[s]))
      
      
    }
    all_means[s] = mean(all_vals)
    all_sds[s] = sd(all_vals)
  }
  
  return(data.frame(winSize = window_sizes, sizeMean=all_means, sizeSD = all_sds))
}


getStatus = function(disease, reference, column='Mut_Prot_Pos', subsetToPLP=T, subsetToPLPVUS=F) {
  ## Given the disease and reference data and a column name,
  ## calculates the affected T/F vector and the one dimensional
  ## locations
  all_pos = sort(unique(c(disease$Mut_Prot_Pos, reference$Mut_Prot_Pos)))
  if (subsetToPLP) {
    aff_pos = unique(subset(disease, Gen_LM_VarClass%in%c("Pathogenic", "Likely Pathogenic"))$Mut_Prot_Pos)
  } else if (subsetToPLPVUS) {
    aff_pos = unique(subset(disease, Gen_LM_VarClass%in%c("Pathogenic", "Likely Pathogenic", "Unknown Significance"))$Mut_Prot_Pos)
  } else {
    aff_pos = unique(disease$Mut_Prot_Pos)
  }
  
  affected = all_pos%in%aff_pos
  
  return(list(affected=affected, position=all_pos))
}

# Gets the Status but ignores any site where a VUS is observed
getStatusExcludeVUS <- function(disease, reference, column = 'Mut_Prot_Pos') {
  all_pos = sort(unique(c(disease$Mut_Prot_Pos, reference$Mut_Prot_Pos)))
  aff_pos = unique(subset(disease, Gen_LM_VarClass%in%c("Pathogenic", "Likely Pathogenic"))$Mut_Prot_Pos)
  vus_pos = unique(subset(disease, Gen_LM_VarClass%in%c("Unknown Significance"))$Mut_Prot_Pos)
  vus_pos <- vus_pos[!vus_pos%in%aff_pos]
  ## Remove vus pos from all pos
  all_pos <- all_pos[!all_pos%in%vus_pos]
  affected <- all_pos %in% aff_pos
  
  return(list(affected = affected, position = all_pos))
}

# Writes the SSS output to a set of files
writeSSSoutput = function(dir, test_stats, perms, prefix="") {
  write.table(file=paste(dir,"/",prefix,"max_stat.txt", sep=""), test_stats$max_stat,quote=F, row.names=F, col.names=F)
  write.table(file=paste(dir,"/",prefix,"test_stats.txt", sep=""), test_stats$test_stats, sep="\t", quote=F, row.names=F, col.names=F)
  write.table(file=paste(dir,"/",prefix,"window_loc.txt", sep=""), test_stats$loc, sep="\n", quote=F, row.names=F, col.names=F)
  write.table(file=paste(dir,"/",prefix,"window_sizes.txt", sep=""), test_stats$sizes, sep="\n", quote=F, row.names=F, col.names=F)
  write.table(file=paste(dir,"/",prefix,"perms.txt", sep=""), perms, sep="\n", quote=F, row.names=F, col.names=F)
}



surface_windows <- function(affected, amino_acid, distance_frame, window_size, target_aa) {
  ## Affected: boolean of the 1/0 status of the amino acids
  ## amino acids: list of the amino acid residues corresponding to each affected value
  ## distance frame: data frame of the distances between amino acids, values pre-calculated
  ## window_size: max. allowable distance between amino acids to be included in window
  ## target_aa: amino acid to target (ie. window center)
  
  
  # Goal: get list of AAs within distance and outside of distance
  in_dist <- subset(distance_frame, (Distance_Strict < window_size)& (AA1==target_aa | AA2 == target_aa))
  aa_in_dist <- unique(c(in_dist$AA1, in_dist$AA2))
  
  aff_in <- affected[amino_acid%in%aa_in_dist]
  
  thislik <- lik_sss(sum(aff_in), sum(affected), length(aff_in), length(affected))
  
  return(thislik)
}

surface_windows_nearest <- function(affected, amino_acid, distance_frame, window_size, target_aa) {
  ## Affected: boolean of the 1/0 status of the amino acids
  ## amino acids: list of the amino acid residues corresponding to each affected value
  ## distance frame: data frame of the distances between amino acids, values pre-calculated - THIS FXN EXPECTS PRE-ORDERED values
  ## window_size: fixed number of closest amino acids to check
  ## target_aa: amino acid to target (ie. window center)
  
  
  # Goal: get list of AAs within distance and outside of distance
  this_aa <- distance_frame[distance_frame$AA1==target_aa | distance_frame$AA2 == target_aa,]
  in_dist <- this_aa[1:(window_size+1),]
  aa_in_dist <- unique(c(in_dist$AA1, in_dist$AA2))
  
  aff_in <- affected[amino_acid%in%aa_in_dist]
  
  thislik <- lik_sss(sum(aff_in), sum(affected), length(aff_in), length(affected))
  
  return(thislik)
}


oneWinSizeSurface = function(affected, amino_acid, distance_frame, window_size, target_aas) {
  # Affected - T/F of affected status
  # Position - Matrix of positions
  # Window locations - List of the window locations
  # Window function - function to calcualte likelihood in windows, default is euclidean distance
  # Window size - size of window
  
  
  resout <- sapply(target_aas, surface_windows, affected=affected, amino_acid=amino_acid, window_size=window_size, distance_frame=subset(distance_frame, Distance_Strict < window_size ))
  
  ## Return the vector of the window test statistics
  return(resout)
}

oneWinSizeSurface_nearest = function(affected, amino_acid, distance_frame, window_size, target_aas) {
  # Affected - T/F of affected status
  # Position - Matrix of positions
  # Window locations - List of the window locations
  # Window function - function to calcualte likelihood in windows, default is euclidean distance
  # Window size - size of window
  # Distance frame - distances between AAs - EXPECTS THIS ORDERERD
  
  resout <- sapply(target_aas, surface_windows_nearest, affected=affected, amino_acid=amino_acid, window_size=window_size, distance_frame=distance_frame)
  
  ## Return the vector of the window test statistics
  return(resout)
}

getMaxLikSurface <- function(affected, amino_acid, distance_frame, window_sizes, target_aas) {
  ## affected - vector of 1/0 affected status for amino acids
  ## amino acid - list of amino acid residues of the corresponding entry in the affected vector
  ## distance_frame - dataframe showing the distances between amino acids, here assuming surface distance
  ## needs columns AA1, AA2, and Distance_Strict
  ## window sizes - vector of window sizes to check
  ## target_aas - amino acids to use as window centers
  allTestStats = matrix(NA, nrow=length(target_aas), ncol=length(window_sizes))
  
  for (s in 1:length(window_sizes)) {
    
    thislik <- oneWinSizeSurface(affected, amino_acid, distance_frame, window_sizes[s], target_aas)
    
    allTestStats[,s] = thislik
  }
  
  return( list(max_stat = max(allTestStats), test_stats = allTestStats, loc = target_aas, sizes = window_sizes))
}

getMaxLikSurface_nearest <- function(affected, amino_acid, distance_frame, window_sizes, target_aas) {
  ## affected - vector of 1/0 affected status for amino acids
  ## amino acid - list of amino acid residues of the corresponding entry in the affected vector
  ## distance_frame - dataframe showing the distances between amino acids, here assuming surface distance
  ## needs columns AA1, AA2, and Distance_Strict
  ## window sizes - vector of window sizes to check (this will be number of closest residues to assess)
  ## target_aas - amino acids to use as window centers
  allTestStats = matrix(NA, nrow=length(target_aas), ncol=length(window_sizes))
  distance_frame <- unique(distance_frame[order(distance_frame$Distance_Strict),])
  for (s in 1:length(window_sizes)) {
    print(s)
    thislik <- oneWinSizeSurface_nearest(affected, amino_acid, distance_frame, window_sizes[s], target_aas)
    
    allTestStats[,s] = thislik
  }
  
  return( list(max_stat = max(allTestStats), test_stats = allTestStats, loc = target_aas, sizes = window_sizes))
}


getAAsInWindow <- function(target_aa, distance_frame, window_size) {
  in_win <- subset(distance_frame, Distance_Strict < window_size & (AA1==target_aa | AA2 == target_aa))
  
  return(unique(c(in_win$AA1, in_win$AA2)))
}

getAAsInWindow_nearest <- function(target_aa, distance_frame, window_size) {
  distance_frame <- unique(distance_frame[order(distance_frame$Distance_Strict),])
  in_win <- distance_frame[1:(window_size+1),]
  
  return(unique(c(in_win$AA1, in_win$AA2)))
}

maxLikSurfacePerms <- function(affected, amino_acid, distance_frame, window_sizes, target_aas, nperms=10) {
  ## Replicates prmutations to get the mpirical distribution of the maximum surface analysis statistic
  r = proc.time()
  
  max_stats <- replicate(nperms, singleSurfacePerm( affected=affected, amino_acid=amino_acid, distance_frame=distance_frame, window_sizes=window_sizes, target_aas = target_aas ))
  
  rtime = proc.time() - r
  print(rtime)
  return(max_stats)
}

singleSurfacePerm <- function(affected, amino_acid, distance_frame, window_sizes, target_aas) {
  thisperm = getMaxLikSurface(sample(affected, length(affected), replace=F), amino_acid, distance_frame, window_sizes, target_aas)
  return(thisperm$max_stat)
}

### Going to add a few functions for the S2 that enforce the symmetry

getPermLiksRep_S2_both = function(affected, position, window_locations, window_sizes, chain, nperms =100, window_func = eucWinLik_sparse) {
  ## Better optimized function for getting the likelihoods from the perumtations, still takes a bit of time though
  ## Inputs follow the inputs of getMaxLik
  ## Needs to enforce the symmetry
  ## Assumes input is sorted by chain then by amino acid residue!!!!!!! Very important to remember this!!!
  print("Did you remember to sort the affected data correctly? Which means first by Chain then by ResidueNo")
  r = proc.time()
  rliks = replicate(nperms, singlePermLiks_S2_both(affected, position, window_locations, window_sizes, chain, window_func =  window_func))
  rtime = proc.time() - r
  print(rtime)
  return(rliks) 
}

singlePermLiks_S2_both = function(affected, position, window_locations, window_sizes, chain, window_func = eucWinLik) {
  ## Helper function for the replicate function above
  ## Tried apply but it doesn't really spee dup the loop (only 4 items or so)
  ## Need to enforce the symmetry here
  ## Assumes input is sorted by chain then by amino acid residue!!!!!!! Very important to remember this!!!
  this_lik = numeric(length(window_sizes))
  aff_a = affected[chain=="A"]
  for (i in 1:length(window_sizes)) {
    # Sample for chainA
    this_samp <- sample(aff_a, length(aff_a))
    this_lik[i] = max(oneWinSize_sparse(c(this_samp, this_samp), position, window_locations, window_size = window_sizes[i], window_func = window_func))
    
  }
  return(max(this_lik))
}
