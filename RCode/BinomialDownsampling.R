#before anything else, you need to setup paths:
#for example
#source("C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/paths.R")
#or, in colab:
#source("GRNP_2020/RCode/pathsGoogleColab.R")

source(paste0(sourcePath, "ButterflyHelpers.R"))


library(stringr)
library(dplyr)
library(qdapTools)

#t is a vector of downsampling values between 0 and 1
binomialDownsampling = function(bug, fractionsToKeep) {
  collapsed = bug %>% group_by(gene) %>% do(countslist=c(.$count))#if you get an error here, you probably defined a variable called "c"...
  
  numGenes = dim(collapsed)[1]
  
  #create a matrix with genes as rows and fraction to keep values as columns
  estTotCounts = matrix(data=0,nrow=numGenes,ncol=length(fractionsToKeep))
  
  print(paste0("Genes: ",numGenes))
  
  
  for (i in 1:numGenes) {
    h = hist(collapsed[[2]][[i]], breaks=seq(0.5, max(collapsed[[2]][[i]])+0.5, by=1), plot = F)$counts
    lh = length(h)
    for (j in 1:length(fractionsToKeep)) {
      ftk = fractionsToKeep[j]
      hd = rep(0,lh)
      for (k in 1:lh) {
        dens = dbinom(1:lh, k, ftk)
        hd = hd + dens*h[k]
      }
      estTotCounts[i,j] = sum(hd)
    }
    
    if ((i %% 1000) == 0) {
      print(i)
    }
  }
  
  #annoying conversion, can probably be done smarter
  nms = c("gene",paste0("d",fractionsToKeep))
  colnames(estTotCounts) = nms[2:length(nms)]
  res = bind_cols(tibble(gene=collapsed$gene), as_tibble(estTotCounts))
  return(res)

}
