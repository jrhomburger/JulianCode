### Admixture Plotting Function ###

## Inputs
# Famfile: Tells what order in the admixture output the individuals are in (can be PLINK famfile format)
# the popinfo file must have columns SampleID, Population, and IndOrder (can change name of 'Population' column in parameters)
# These will be used for the plotting and ordering of the data
# The Admixtable is the Q File output from ADMIXTURE
plotAdmixture = function(famfile, popinfo, admixtable, colors=c("dodgerblue4", "darkgreen", "firebrick", "orange", "yellowgreen", "darkolivegreen4", "steelblue"), poplabels="Population", cex.names=0.7 ,las=2, adj=0.5) {
  famfile$index = 1:nrow(famfile)
  popinfo = subset(popinfo, SampleID%in%famfile$V2)
  fampop = merge(famfile, popinfo, by.x="V2", by.y="SampleID")
  fampop = fampop[order(fampop$index),]
  k = ncol(admixtable)
  if (k > length(unique(colors))) {
    print("Warning: more clusters than specified colors.")
  }
  
  knames = paste("K", 1:k, sep="")
  kcombo = cbind(admixtable, fampop[[poplabels]], fampop$IndOrder, fampop$V2)
  names(kcombo) = c(knames, "Population", "IndOrder", "SampleID")
  kcombo = kcombo[order(kcombo$IndOrder),]
  # Make the combine file
  
  spaces <- c(0,diff(kcombo$Population))
  spaces[spaces!=0] <- 2
  index <- c(1, diff(kcombo$Population))
  index[index !=0] <- 1
  names <- kcombo$Population
  borders <- seq(1, length(index))[index==1]
  round(diff(c(borders,length(index) ) )/2)
  offset <- round(diff(c(borders,length(index) ) )/2)
  newnames <- names
  newnames <- rep("", length(names))
  newnames[borders+offset] <- as.character(names[borders+offset])
  
  barplot(t(kcombo[,1:k]), col=colors, border=NA, ylab= paste("K=",k,sep=""), space=spaces, cex.names=cex.names ,las=las, adj=adj, names=newnames)
  
}

