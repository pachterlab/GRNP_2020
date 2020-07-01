#before anything else, you need to setup paths:
#for example
#source("C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/paths.R")
#or, in colab:
#source("GRNP_2020/RCode/pathsGoogleColab.R")

library(magrittr)
#library(data.table)
#library(DropletUtils)
library(Matrix)

library(dplyr)
#library(textTinyR)

#library(Rcpp)

namedSave <- function(objlist, objNames, file) {
  names(objlist) <- objNames
  save(list=names(objlist), file=file, envir=list2env(objlist))
}



genBugFileName = function(dsFactor) {
  paste0("Bug_", dsFactor*100, ".RData")
}

genBugObjName = function(name, dsFactor) {
  paste0("bug_", name, "_", dsFactor*100)
}

getBug = function(name, dsFactor = 1) {
  get(genBugObjName(name, dsFactor), envir=.GlobalEnv)
}

getStats = function(name, dsFactor = 1) {
  get(paste0("stats", name), envir=.GlobalEnv)
}

loadBug <- function(name, dsFactor = 1, fig_data_path = figure_data_path) {
  subpath = paste0(fig_data_path, name, "/")
  load(paste0(subpath, genBugFileName(dsFactor)), envir=.GlobalEnv)
}

rmBug <- function(name, dsFactor = 1) {
  rm(list=genBugObjName(name, dsFactor), envir=.GlobalEnv)
}

loadStats <- function(name, fig_data_path = figure_data_path) {
  subpath = paste0(fig_data_path, name, "/")
  load(paste0(subpath, "Stats.RData"), envir=.GlobalEnv)
}

#down-samples a bunch of times and gets the mean histograms from that
#saves the result to file
#expects that the bug is loaded
generatePooledHistogramDS <- function(dsid, fractionToKeep = 0.1, numSamplings=10) {
  bug = getBug(dsid, 1)
  hs = downSampleManyTimesAndGetHist(bug, fractionToKeep, numSamplings)
  subpath = paste0(figure_data_path, dsid, "/")
  namedSave(list(hs), list(statsName = paste0("pooledHistDS",dsid)), paste0(subpath, "pooledHistDS.RData"))
}

loadPooledHistogramDS <- function(dsid) {
  subpath = paste0(figure_data_path, dsid, "/")
  load(paste0(subpath, "pooledHistDS.RData"), envir=.GlobalEnv)
}

generatePooledHistogram <- function(dsid) {
  bug = getBug(dsid, 1)
  hs = getDsHist(bug)
  subpath = paste0(figure_data_path, dsid, "/")
  namedSave(list(hs), list(statsName = paste0("pooledHist",dsid)), paste0(subpath, "pooledHist.RData"))
}

loadPooledHistogram <- function(dsid) {
  subpath = paste0(figure_data_path, dsid, "/")
  load(paste0(subpath, "pooledHist.RData"), envir=.GlobalEnv)
}

downSampleBUG <- function(bug, fractionToKeep) {
  print("1")
  totCounts = sum(bug$count)
  print("2")
  toSample = totCounts*fractionToKeep
  print("3")
  indices = sample(totCounts, toSample, replace = F)
  print("4")
  borders = c(0.5,cumsum(bug$count) + 0.5)
  print("5")
  newCounts = hist(indices, borders, plot=F)
  print("6")
  sel = newCounts$count != 0
  print("7")
  bugDS = bug[sel,]
  print("8")
  bugDS$count = newCounts$count[sel]
  print("9")
  return(bugDS)
}

#Downsamples numTimes times and sums up all the UMI counts
downSampleBUGNTimes <- function(bug, fractionToKeep, numTimes) {
  print("Preparing...")
  totCounts = sum(bug$count)
  toSample = totCounts*fractionToKeep
  borders = c(0.5,cumsum(bug$count) + 0.5)
  UMIs = matrix(0, nrow=length(unique(bug$gene)), ncol=numTimes)
  firstRound = T
  for (i in 1:numTimes) {
    print(i)
    indices = sample(totCounts, toSample, replace = F)
    newCounts = hist(indices, borders, plot=F)
    #sel = newCounts$count != 0
    #bugDS = bug[sel,] #don't bother fixing the counts per UMI in bugDS, we don't need them
    bugDS = bug
    bugDS$count = newCounts$counts
    tibb = bugDS %>% group_by(gene) %>% summarize(n=length(count[count != 0]))
    UMIs[,i] = tibb$n
    #get counts
    if (firstRound) {
      firstRound = F
      genes = tibb$gene
    }
  }

  sumCounts = rowSums(UMIs)

  #create a tibble from the matrix and the genes
  res = tibble(gene=genes,n=sumCounts)


  return(res)
}

#gets the histogram for each gene
getDsHist <- function(bug) {
  collapsed = bug %>% group_by(gene) %>% do(countslist=c(.$count))
  hgram = matrix(0, nrow=length(collapsed$gene), ncol=100)

  for (j in 1:length(collapsed$gene)) {
    #include the zeros that were not filtered above, but skip them below
    h = hist(collapsed[[j,2]], breaks=c(seq(-0.5, 99.5, by=1), 10000000), plot = F)
    #ignore the number of zeros, we just keep them in to preserve the gene order
    hgram[j,] = h$counts[-1]
  }

  row.names(hgram) = collapsed$gene

  return(hgram)
}

#gives a histogram having the number of copies, 1, 2 and >2
#it is built from down-sampling the same bug many times
#use 100 as histogram limit - all counts above 100 ends up at 100
downSampleManyTimesAndGetHist <- function(bug, fractionToKeep, numTimes=10) {
  print("Preparing...")
  totCounts = sum(bug$count)
  toSample = totCounts*fractionToKeep
  borders = c(0.5,cumsum(bug$count) + 0.5)
  hgram = matrix(0, nrow=length(unique(bug$gene)), ncol=100)
  for (i in 1:numTimes) {
    print(i)
    indices = sample(totCounts, toSample, replace = F)
    print(paste0(i,":1"))
    newCounts = hist(indices, borders, plot=F)
    print(paste0(i,":2"))
    #sel = newCounts$count != 0
    #bugDS = bug[sel,] #don't bother fixing the counts per UMI in bugDS, we don't need them
    bugDS = bug
    bugDS$count = newCounts$counts
    print(paste0(i,":3"))
    collapsed = bugDS %>% group_by(gene) %>% do(countslist=c(.$count))
    print(paste0(i,":4"))
    for (j in 1:length(collapsed$gene)) {
      #include the zeros that were not filtered above, but skip them below
      h = hist(collapsed[[j,2]], breaks=c(seq(-0.5, 99.5, by=1), 10000000), plot = F)
      #ignore the number of zeros, we just keep them in to preserve the gene order
      hgram[j,] = hgram[j,] + h$counts[-1]
    }
    print(paste0(i,":5"))
  }

  row.names(hgram) = collapsed$gene
  #filter out rows with 0 counts
  hgram = hgram[rowSums(hgram) != 0,]

  return(hgram)
}



totalCPUHistogram <- function(bug) {
  return (hist(bug$count, breaks=seq(0.5, max(bug$count)+0.5, by=1), plot = F))
}

#so, t here is not the t in the Good-Toulmin formula, but for how many times we should predict forward
goodToulmin <- function(histo, t) {
  res = 0
  for (i in 1:length(histo$counts)) {
    res = res - ((-(t-1))^i)*histo$counts[i]
  }
  #also add the already known reads:
  res = res + sum(histo$counts)
  return (res)
}



#uses quantile normalization
poolHistograms <- function(dsid, dsBug, otherIds, useDownSampledPools = T) {

  histLength = 100

  collapsed = dsBug %>% group_by(gene) %>% do(countslist=c(.$count))
  umis = dsBug %>% group_by(gene) %>% summarize(UMIs = n())
  hgramDensity = matrix(0, nrow=length(collapsed$gene), ncol=histLength)
  for (i in 1:length(collapsed$gene)) {
    h = hist(collapsed[[i,2]], breaks=c(seq(-0.5, 99.5, by=1), 10000000), plot = F)
    #ignore the number of zeros, we just keep them in to preserve the gene order
    hgramDensity[i,] = h$density[-1]
  }
  row.names(hgramDensity) = collapsed$gene

  #loop through all the other datasets and add them
  numGenes = length(collapsed$gene)#only work with the genes that exist in dsid, the other ones are not of any interest
  numDs = length(otherIds)

  allHistsBefNorm = vector(mode = "list", length = numDs)

  for(d in 1:numDs) {
    if (useDownSampledPools) {
      phOther = get(paste0("pooledHistDS",otherIds[[d]]), envir=.GlobalEnv)
    } else {
      phOther = get(paste0("pooledHist",otherIds[[d]]), envir=.GlobalEnv)
    }
    umisOther = tibble(gene = row.names(phOther), UMIs = rowSums(phOther))

    phOther = phOther / rowSums(phOther)

    #make sure the rows come in the same order in each dataset (i.e. that the genes are matched)
    ind = match(row.names(hgramDensity), row.names(phOther))
    newPhOther = hgramDensity
    newPhOther[!is.na(ind),] = phOther[ind[!is.na(ind)],]
    newPhOther[is.na(ind),] = rep(0,histLength)
    allHistsBefNorm[[d]] = newPhOther

    umis = left_join(umis, umisOther, by="gene")
  }

  #make sure the umis are in the same order as the histograms, i.e. that the genes are matched
  umis = umis[match(umis$gene, row.names(hgramDensity)),]#there should be no NAs here

  normDensityHists = allHistsBefNorm #allocation only

  #Now normalize by quantile normalization
  for(d in 1:numDs) {
    dsHist = normDensityHists[[d]]
    for(histInd in 1:histLength) {
      vals = sort(hgramDensity[,histInd])
      ord = sort(dsHist[,histInd], index.return=T)
      dsHist[ord$ix,histInd] = vals
    }
    normDensityHists[[d]] = dsHist
  }

  return (list(umis=umis, hists = c(list(hgramDensity),normDensityHists)))
}

