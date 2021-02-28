#before anything else, you need to setup paths:
#for example
#source("C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/paths.R")
#or, in colab:
#source("GRNP_2020/RCode/pathsGoogleColab.R")

library(preseqR)
library(magrittr)
library(Matrix)

source(paste0(sourcePath,"modZTNB.R"))


#t can be a vector, it is the prediction range (i.e. 2 means predict to the double amount of reads)
upSampleAndGetMeanExprPreSeqZTNB <- function(bugFile, t, incTol = 1e-5, iterIncTol = 200) {
  
  collapsed = bugFile %>% group_by(gene) %>% do(countslist=c(.$count))#if you get an error here, you probably defined a variable called "c"...
  
  numGenes = dim(collapsed)[1]
  
  estTotCounts = matrix(data=0,nrow=numGenes,ncol=length(t))
  
  print(paste0("Genes: ",numGenes))
  
  
  for (i in 1:numGenes) {
    h = hist(collapsed[[2]][[i]], breaks=seq(0.5, max(collapsed[[2]][[i]])+0.5, by=1), plot = F)
    freq = h$mids
    counts = h$counts
    added = 0
    #preseq cannot handle if we have only ones, so modify the histogram slightly
    if ((length(freq)==1) & (freq[1] == 1)) {
      added = 2
      freq = c(1,2)
      counts = c(counts[1]+1,1)#room for improvement here
    }
    dd = as.matrix(data.frame(freq,counts));
    rSAC = mod.ztnb.rSAC(dd, incTol = incTol, iterIncTol = iterIncTol);
    newCounts = rSAC(t)
    newCounts[newCounts < 0] = 0
    if ((i %% 1000) == 0) {
      print(i)
    }
    estTotCounts[i,] = newCounts - added;
  }
  
  #annoying conversion, can probably be done smarter
  nms = c("gene",paste0("p",t))
  colnames(estTotCounts) = nms[2:length(nms)]
  res = bind_cols(tibble(gene=collapsed$gene), as_tibble(estTotCounts))
  return(res)
}

#t can be a vector, it is the prediction range (i.e. 2 means predict to the double amount of reads)
upSampleAndGetMeanExprPreSeqDs <- function(bugFile, t, mt=2) {
  
  collapsed = bugFile %>% group_by(gene) %>% do(countslist=c(.$count))#if you get an error here, you probably defined a variable called "c"...
  
  numGenes = dim(collapsed)[1]
  
  estTotCounts = matrix(data=0,nrow=numGenes,ncol=length(t))
  
  print(paste0("Genes: ",numGenes))
  
  
  for (i in 1:numGenes) {
    h = hist(collapsed[[2]][[i]], breaks=seq(0.5, max(collapsed[[2]][[i]])+0.5, by=1), plot = F)
    freq = h$mids
    counts = h$counts
    added = 0
    #preseq cannot handle if we have only ones, so modify the histogram slightly
    if ((length(freq)==1) & (freq[1] == 1)) {
      added = 2
      freq = c(1,2)
      counts = c(counts[1]+1,1)#room for improvement here
    }
    dd = as.matrix(data.frame(freq,counts));
    rSAC = ds.rSAC(dd, mt=mt)
    newCounts = rSAC(t)
    newCounts[newCounts < 0] = 0
    if ((i %% 1000) == 0) {
      print(i)
    }
    estTotCounts[i,] = newCounts - added;
  }
  
  #annoying conversion, can probably be done smarter
  nms = c("gene",paste0("p",t))
  colnames(estTotCounts) = nms[2:length(nms)]
  res = bind_cols(tibble(gene=collapsed$gene), as_tibble(estTotCounts))
  return(res)
}


#predicts filling in lowly expressed genes from other datasets
poolPrediction <- function(bugFile, t=10, poolHistList, usePoolLimit = 100) {
  mt=2
  grouped = bugFile %>% group_by(gene)
  collapsed = grouped %>% do(countslist=c(.$count))#if you get an error here, you probably defined a variable called "c"...
  umis = grouped %>% summarize(umis=n())

  UMIs = poolHistList[[1]] #the UMIs from the dataset to predict are included here
  hists = poolHistList[[2]]

  numDs = length(hists)

  #now, use the UMIs to calculate how much weight that should be given to each dataset for each gene:
  UMIs[is.na(UMIs)] = 0
  totUMIs = rowSums(UMIs[,-1])

  weights = UMIs[,-1]/totUMIs
  #rowSums(weights)#check, should all be 1, ok

  #Now, multiply each dataset with its weight for each gene
  weighted = hists
  for (i in 1:numDs) {
    weighted[[i]] = weighted[[i]] * weights[,i] # so, R will multiply each row with the same value in the vector, which is what we want
  }

  #now we add the weighted histograms to get a weighted mean
  totHist = weighted[[1]]
  for(i in 2:numDs) {
    totHist = totHist + weighted[[i]]
  }

  #since they are quantile normalized, the density may not add up to exactly 1, so scale that
  totHistDensity = totHist / rowSums(totHist)

  totHist = totHistDensity * UMIs[[2]]

  #loop through all the other datasets and add them
  numGenes = length(UMIs$gene)#only work with the genes that exist in dsid, the other ones are not of any interest
  estTotCounts = matrix(data=0,nrow=numGenes,ncol=length(t))
  print(paste0("Genes: ",numGenes))

  for (i in 1:numGenes) {
    if (umis[i,2] < usePoolLimit) {
      counts = totHist[i,]
    } else {
      counts = hists[[1]][i,]*UMIs[i,2, drop=T]
    }

    #remove trailing zeros
    for (x in 100:1) {
      if (counts[x] != 0) {
        break;
      }
    }
    counts = counts[1:x]
    freq = 1:x

    #preseq cannot handle if we have only ones, so modify the histogram slightly
    added = 0
    if (counts[1] == 0) {
      added = added + 1
      counts[1] = 1
    }
    if (length(counts) == 1) {
      freq = c(1,2)
      counts = c(counts[[1]], 1)
      added = added + 1
    }
    if (counts[2] == 0) {
      added = added + 1
      counts[2] = 1
    }

    dd = as.matrix(data.frame(freq,counts));
    rSAC = ds.rSAC(dd, mt=mt)
    #rSAC = mod.ztnb.rSAC(dd, incTol = 1e-5, iterIncTol = 200);
    newCounts = rSAC(t)
    newCounts[newCounts < 0] = 0
    if ((i %% 1000) == 0) {
      print(i)
    }
    estTotCounts[i,] = newCounts - added
  }

  #annoying conversion, can probably be done smarter
  nms = c("gene",paste0("p",t))
  colnames(estTotCounts) = nms[2:length(nms)]
  res = bind_cols(tibble(gene=collapsed$gene), as_tibble(estTotCounts))
  return(res)
}

#t is a multiplier of the counts, i.e. 2 means multiplying to the double
predPreSeqDS <- function(histo, t, mt) {
  freq = histo$mids
  counts = histo$counts
  #if only ones, add one extra umi with two copies, the algorithm cannot handle it otherwise
  if (length(counts) == 1 & freq[1] == 1) {
    counts = c(counts,1)
    freq = c(1,2)
  }

  dd = as.matrix(data.frame(freq,counts));
  rSAC = ds.rSAC(dd, mt=mt)
  return (rSAC(t))
}

#This has been modified to use the right function that is reasonably fast
#t is a multiplier of the counts, i.e. 2 means predicting at the double number of counts
predPreSeqZTNB <- function(histo, t, incTol = 1e-5, iterIncTol = 200) {
  freq = histo$mids
  counts = histo$counts
  #if only ones, add one extra umi with two copies, the algorithm cannot handle it otherwise
  if (length(counts) == 1 & freq[1] == 1) {
    counts = c(counts,1)
    freq = c(1,2)
  }
  dd = as.matrix(data.frame(freq,counts));
  rSAC = mod.ztnb.rSAC(dd, incTol = incTol, iterIncTol = iterIncTol)
  return (rSAC(t))
}


#t is a multiplier of the counts, i.e. 2 means predicting at the double number of counts
#this is the "best practice" function
predPreSeq <- function(histo, t, mt) {
  #h = hist(bugFile$V3, breaks=seq(min(bugFile$V3)-0.5, max(bugFile$V3)+0.5, by=1), xlim=c(0,30), xlab="Copies per UMI", plot = F)
  freq = histo$mids
  counts = histo$counts
  #if only ones, add one extra umi with two copies, the algorithm cannot handle it otherwise
  if (length(counts) == 1 & freq[1] == 1) {
    counts = c(counts,1)
    freq = c(1,2)
  }
  dd = as.matrix(data.frame(freq,counts));
  rSAC = mod.preseqR.rSAC.fixed(dd, mt=mt)#preseqR.rSAC(dd, mt=mt)
  return (rSAC(t))
}


