library(preseqR)
library(BUSpaRse)
library(TENxBUSData)
library(ggplot2)
library(magrittr)
library(data.table)
library(Seurat)
library(DropletUtils)
library(Matrix)

library(dplyr)
library(textTinyR)
#install.packages("DescTools")
library(DescTools)


sourcePath = "GRNP_2020/RCode/"
dataPath = "data/"
figure_data_path = "figureData/"
figure_path = "figures/"
source(paste0(sourcePath,"modZTNB.R"))



createTr2g <- function(fasta_file, kallisto_out_path) {
  tr2g <- transcript2gene(fasta_file = fasta_file,
                               kallisto_out_path = kallisto_out_path)
  write.table(tr2g, paste0(kallisto_out_path, "/transcripts_to_genes.txt"), quote=F, row.names = F, col.names=F, sep="\t")
}

#gives you the standard ggplot colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#used for investigating single-copy molecules
ClosestDists <- function(bug, subBug, UMILength) {
  #make sure subbug is max 2000 UMIs to shorten the execution time:
  if (dim(subBug)[1] > 2000) {
    print("Down-sampling to 2000 UMIs")
    samp = sort(sample(1:(dim(subBug)[1]), 2000))
    subBug = subBug[samp,]
  }

  currCell = ""
  dists = rep(0,UMILength)
  print(paste0("Will process ", dim(subBug)[1], " UMIs"))
  for (i in 1:(dim(subBug)[1])) {
    if (i%%100 == 0) {
      print(i)
    }
    if (subBug$barcode[i] != currCell) {
      currCell = subBug$barcode[i];
      currCellBug = bug[bug$barcode==currCell,]
    }
    currUMI = subBug$UMI[[i]]
    bestDist = UMILength;
    for (j in 1:(dim(currCellBug)[1])) {
      dist = stringdist::stringdist(currCellBug$UMI[[j]], currUMI, method = "hamming")
      if (dist > 0 & dist < bestDist) {
        bestDist = dist;
        if (bestDist == 1){
          break
        }
      }
    }
    dists[bestDist] = dists[bestDist] + 1;
  }
  return (dists)
}

#convert gene indices to gene symbols
library(qdapTools)
geneIndices2Symbols <- function(geneIndices, genesFile, Tr2gFile) {
  genenames = read.table(genesFile, stringsAsFactors = F)$V1
  tr2g = read.table(Tr2gFile, stringsAsFactors = F)
  lookupTable = tr2g[,2:3]
  lookupTable= unique(lookupTable)
#  lookupTableTranscr = tr2g[,c(1,3)]
  numGenes = length(geneIndices)

  inGenes = geneIndices;#allocate

  for (i in 1:numGenes) {
    inGenes[i] = genenames[1 + as.numeric(geneIndices[i])]
  }
  #check that the ens genes are really unique in the lookup table:
#  length(lookupTable$V2)
#  length(unique(lookupTable$V2)) #ok

  outGenes = lookup(inGenes, lookupTable)

  return (outGenes)
}

library(stringr)

#dir should include slash at the end
readBug <- function(dir) {
  print(paste0("Reading BUG from ", dir, " ..."))
  bug = read.table(paste0(dir,"bus_output/bug.txt"), stringsAsFactors = F)
  #gllengths = sapply(ngbug[,3],function(s) str_count(s, fixed(",")) + 1)
  print("Filtering multi-mapped reads...")
  isMult = sapply(bug[,3],function(s) grepl(",",s, fixed=T))

  print (paste0("Fraction multi-mapped reads: ", sum(isMult) / dim(bug)[1]))

  uniquelymapped = bug[!isMult,] #just ignore the multimapped reads

  print("Converting genes...")

  #now convert gene index to gene id
  uniquelymapped$V3 = geneIndices2Symbols(uniquelymapped$V3, paste0(dir,"bus_output/coll.genes.txt"), paste0(dir,"bus_output/transcripts_to_genes.txt"))

  colnames(uniquelymapped) = c("barcode", "UMI", "gene", "count")

  print("Done")

  return (uniquelymapped)
}

downsampleBugs <- function(bug, fracs) {
  print(paste0("Down-sampling in total ", length(fracs), " bugs:"))
  bugs = vector(mode = "list", length = length(fracs))
  for (i in 1:length(fracs)) {
    print(paste0(i, ": Down-sampling to ", fracs[[i]]))
    if (fracs[[i]] == 1) {
      bugs[[i]] = bug
    } else {
      bugs[[i]] = downSampleBUG(bug,fracs[[i]])
    }
  }
  print("Done")
  return(bugs)
}

formatFracName<- function(name, frac, type) {
  paste0(type, "_", name,"_d_", frac*100)
}

fracOnesFunc <- function(d) {
  sum(d == 1) / length(d)
}


getStatsFromBugs <- function(dsBugs, fracs, name) {
  stats = vector(mode = "list", length = length(fracs))
  for (i in 1:length(dsBugs)) {
    tmp = dsBugs[[i]] %>% group_by(gene) %>% summarise(UMIs = n(), counts=sum(count), CPM = n(), fracOnes = fracOnesFunc(count), countsPerUMI = mean(count))
    #CPM needs to be fixed, only UMI counts right now
    tmp$CPM = tmp$CPM*10^6/sum(tmp$CPM)
    #set the right column names here
    cn = c("gene",
           formatFracName(name, fracs[[i]], "UMIs"),
           formatFracName(name, fracs[[i]], "Counts"),
           formatFracName(name, fracs[[i]], "CPM"),
           formatFracName(name, fracs[[i]], "FracOnes"),
           formatFracName(name, fracs[[i]], "CountsPerUMI")
    )
    colnames(tmp) = cn
    stats[[i]] = tmp
  }

  #now merge to one table
  stat = stats[[1]]
  if (length(fracs) > 1) {
    for(i in 2:length(fracs)) {
      stat = inner_join(stat, stats[[i]], by="gene")
    }
  }

  return (stat)
}


namedSave <- function(objlist, objNames, file) {
  names(objlist) <- objNames
  save(list=names(objlist), file=file, envir=list2env(objlist))
}


createStandardBugsData <- function(bugdir, name, fracs, UmisPerCellLimit = 200, fig_data_path = figure_data_path) {
  #Generate data
  print(paste0("Generating data for ", name))
  bug = readBug(bugdir)
  #filter out low quality cells
  #Should have more than 200 UMIs
  UMIsPerCell = bug %>% group_by(barcode) %>% tally()
  sum(UMIsPerCell$n > UmisPerCellLimit)
  filtBug = bug[bug$barcode %in% UMIsPerCell$barcode[UMIsPerCell$n > UmisPerCellLimit],]

  #skip mitochondrial content for now; don't think it matters for this application
  #should have mitochondrial content of less than 10%
  #genelist = filtBug %>% group_by(barcode) %>% do(genes=c(.$gene))

  #mc = rep(0,length(genelist$barcode))

  #for (i in 1:length(genelist)) {
  #  tmp <- grep(pattern = "^mt-", x = genelist$genes[[i]], value = TRUE)
  #  mc[i] = length(tmp) /length(genelist$genes[[i]])
  #}


  dsBugs = downsampleBugs(filtBug, fracs)
  stats = getStatsFromBugs(dsBugs, fracs, name)

  #create folder
  subpath = paste0(fig_data_path, name, "/")
  if (!file.exists(paste0(fig_data_path, name))){ #there must be no slash at the end here...
    dir.create(file.path(fig_data_path, name))
  }

  #generate proper variable names before saving
  fracsName = paste0("fracs",name)
  bugName = paste0("bug",name)
  bugsName = paste0("dsBugs",name)
  statsName = paste0("stats",name)

  print("Saving BUG...")
  namedSave(list(filtBug, fracs), list(bugName, fracsName), paste0(subpath, "Bug.RData"))
#  save(get(bugName), get(fracsName), file=paste0(subpath, "Bug.RData"))
  print("Saving BUGs...")
  namedSave(list(dsBugs, fracs), list(bugsName, fracsName), paste0(subpath, "DsBugs.RData"))
  #save(get(bugsName), get(fracsName), file=paste0(subpath, "DsBugs.RData"))
  print("Saving stats...")
  namedSave(list(stats), list(statsName), paste0(subpath, "Stats.RData"))
  #save(list, file=paste0(subpath, "Counts.RData"))
  #save(get(cpmsName), file=paste0(subpath, "Cpms.RData"))
}


loadBug <- function(name, fig_data_path = figure_data_path) {
  subpath = paste0(fig_data_path, name, "/")
  load(paste0(subpath, "Bug.RData"), envir=.GlobalEnv)
}

loadBugs <- function(name, fig_data_path = figure_data_path) {
  subpath = paste0(fig_data_path, name, "/")
  load(paste0(subpath, "DsBugs.RData"), envir=.GlobalEnv)
}

loadStats <- function(name, fig_data_path = figure_data_path) {
  subpath = paste0(fig_data_path, name, "/")
  load(paste0(subpath, "Stats.RData"), envir=.GlobalEnv)
}

#down-samples a bunch of times and gets the mean histograms from that
#saves the result to file
#expects that the bug is loaded
generatePooledHistogramDS <- function(dsid, fractionToKeep = 0.1, numSamplings=10) {
  bug = get(paste0("bug",dsid), envir=.GlobalEnv)
  hs = downSampleManyTimesAndGetHist(bug, fractionToKeep, numSamplings)
  subpath = paste0(figure_data_path, dsid, "/")
  namedSave(list(hs), list(statsName = paste0("pooledHistDS",dsid)), paste0(subpath, "pooledHistDS.RData"))
}

loadPooledHistogramDS <- function(dsid) {
  subpath = paste0(figure_data_path, dsid, "/")
  load(paste0(subpath, "pooledHistDS.RData"), envir=.GlobalEnv)
}

generatePooledHistogram <- function(dsid) {
  bug = get(paste0("bug",dsid), envir=.GlobalEnv)
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


#t can be a vector, it is the prediction range (i.e. 2 means predict to the double amount of reads)
upSampleAndGetMeanExprPreSeqZTNB <- function(bugFile, t, incTol = 1e-5, iterIncTol = 200) {
  
  collapsed = bugFile %>% group_by(gene) %>% do(countslist=c(.$count))#if you get an error here, you probably defined a variable called "c"...
  
  numGenes = dim(collapsed)[1]
  
  estTotCounts = matrix(data=0,nrow=numGenes,ncol=length(t))
  
  print(paste0("Genes: ",numGenes))
  
  
  for (i in 1:numGenes) {
    h = hist(collapsed[[i,2]], breaks=seq(0.5, max(collapsed[[i,2]])+0.5, by=1), plot = F)
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



library(extraDistr)


getCCC = function(x, y) {
  CCC(x,y)$rho.c$est
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

