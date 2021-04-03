#
# 
#

#before anything else, you need to setup paths:
#for example
#source("C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/paths.R")
#or, in colab:
#source("GRNP_2020/RCode/pathsGoogleColab.R")

source(paste0(sourcePath,"ButterflyHelpers.R"))
source(paste0(sourcePath,"ggplotHelpers.R"))
source(paste0(sourcePath,"preseqHelpers.R"))

library(ggplot2)
library(ggpubr)
#library(DESeq2")
library(Seurat)


#normal values, when not varying the other factors
#depth = 1 #(default in function, no need to specify that)
#standard mu is to use the mu distribution mus
#deFC = 4 #(default in function, no need to specify that)
standardExpr = 1000;



#We generate data as negative binomial

#We need to find realistic parameterizations for the negative binomials
#We use a size of 1, which is a reasonable number when looking at highly expressed genes
#The relation between mu and FSCM is then calculated. For each gene, we randomly select
#a fraction of single-copy molecules (FSCM) from the histogram of PBMC_V3_3 (as shown in Fig S1).
#That value is translated to an NB mean value using linear interpolation

calcFSCM = function(mu, size=1) {
  dnbinom(1, mu=mu, size=size) / (1-dnbinom(0, mu=mu, size=size))
}


means = 2^seq(-7,3,by=0.2)
FSCMs = rep(NA, length(means))
for (i in 1:length(means)) {
  FSCMs[i] = calcFSCM(means[i])
}

#Get a list of FSCM values from the PBMC_V3_3 dataset, and convert them to mu
dsid = "PBMC_V3_3"
loadStats(dsid)
statsV3 = getStats(dsid)
FSCMsV3 = statsV3$FracOnes_PBMC_V3_3_d_100[statsV3$UMIs_PBMC_V3_3_d_100 >= 30]
CU = statsV3$CountsPerUMI_PBMC_V3_3_d_100[statsV3$UMIs_PBMC_V3_3_d_100 >= 30]
max(CU) # 8.840637 it makes sense to assume that mean is not larger than 8
mus = approx(x = FSCMs, y = means, xout = FSCMsV3, method = "linear", yleft=min(means), yright=max(means))$y


#This is how the linear interpolation is done
#newMeans = approx(x = FSCMs, y = means, xout = c(0.1, 0.5, 0.8), method = "linear")
#9.018959 1.000000 0.250000
#Test 1
#calcFSCM(0.25) #0.8, ok

#For each gene expression level, we generate pairs of genes with randomly selected amplification + some noise
#A fraction of those are DE. We then make a DE test with and without correction


#returns a histogram and a vector of UMIs per cell
genGeneData = function(numCells, geneExpr, mu, totMoleculesPerCell=12000) {
  h = rep(NA, 1000) #histogram, skip the zeros, i.e. starts at 1
  
  totGeneMolPerCell = geneExpr/10^6*totMoleculesPerCell
  numUMIs = round(totGeneMolPerCell*numCells)
  #Generate CU histogram
  molecules = 0
  #don't accept zero molecules, it messes everything up in the analysis further down
  while (sum(molecules) == 0) {
    molecules = rnbinom(numUMIs, mu=mu, size=1)
  }
  molecules[molecules > 1000] = 1000 #remove outliers for practical reasons
  seenMolecules = molecules[molecules > 0]
  h = hist(seenMolecules, breaks=seq(0.5, 1000+0.5, by=1), plot = F)$counts

  nonZeroMol = sum(h)
  
  #Generate counts across cells distribution
  x = sample(numCells,nonZeroMol,replace=TRUE)
  m = hist(x, breaks=seq(0.5, numCells+0.5, by=1), plot = F)$counts


  return(list(m,h))
}

#some params for debugging
numCells = 2000
fractionDE = 0.5
numGenes = 500

DEExprFactor = 4
geneExpr = 1000
fixedMu=0.5
doCorrect = TRUE

trimZeros = function(h) {
  return(h[1:max( which( h != 0 ))])  
}
#Test 2
#trimZeros(c(1,0,2,3,0,1,0,0,0,0))#1 0 2 3 0 1, ok

predZTNB = function(h, t) {
  #Get the histogram
  freq = 1:length(h)
  counts = h
  added = 0
  #preseq cannot handle if we have only ones, so modify the histogram slightly
  if ((length(freq)==1) & (freq[1] == 1)) {
    added = 2
    freq = c(1,2)
    counts = c(counts[1]+1,1)#room for improvement here
  }
  dd = as.matrix(data.frame(freq,counts));
  rSACZTNB = mod.ztnb.rSAC(dd, incTol = 1e-5, iterIncTol = 200);
  newCountsZTNB = rSACZTNB(t)
  newCountsZTNB[newCountsZTNB < 0] = 0
  return(newCountsZTNB - added);
}

binomialDowns = function(h, fractionToKeep) {
  
  lh = length(h)
  hd = rep(0,lh)
  for (k in 1:lh) {
    dens = dbinom(1:lh, k, fractionToKeep)
    hd = hd + dens*h[k]
  }
  estTotCounts = sum(hd)
}



#mus can either be a single value or a vector, from which one randomly will be chosen
#set mus2 to NA if we want to use the same mu for both conditions (can be used for depth calculations for example)
evaluateCondition = function(geneExpr, mus1, mus2, depth1 = 1, depth2 = 1, DEExprFactor = 4, numCells = 5000, numGenes = 500, fractionDE = 0.5, useBinDowns = FALSE, minLFC = 0.25) {
  geneExpr2 = geneExpr * DEExprFactor

  print("Generating data...")
  
  #randomize which genes should be DE - place the first half in the first group, the other half in the other, no overlaps
  numDEDiv2 = round(numGenes*fractionDE/2) #half increases first group, half the second
  numDE = numDEDiv2*2
  deGenes = sample(numGenes, numDE, replace=FALSE)
  isDE1 = 1:numGenes %in% deGenes[1:numDEDiv2]
  isDE2 = 1:numGenes %in% deGenes[(1:numDEDiv2)+numDEDiv2]
  isDE = isDE1 | isDE2
  #sum(isDE) #ok
  ge1 = rep(geneExpr, numGenes)
  ge1[isDE1] = geneExpr2
  ge2 = rep(geneExpr, numGenes)
  ge2[isDE2] = geneExpr2
  
  #the first numCells cells is the first group, the next numCells is the other group
  geMatrix = matrix(nrow=numGenes, ncol=numCells*2)
  colnames(geMatrix) = as.character(1:(2*numCells))
  rownames(geMatrix) = as.character(1:(numGenes))
  
  histGroup1 = matrix(nrow=numGenes, ncol=1000)
  histGroup2 = matrix(nrow=numGenes, ncol=1000)
  
  meta = data.frame(group=as.factor(c(rep(0,numCells), rep(1,numCells))))
  rownames(meta) = as.character(1:(2*numCells))
  
  mus1Rnd = sample(length(mus1), numGenes, replace=TRUE)
  musr1 = mus1[mus1Rnd]*depth1
  if (is.na(mus2[1])) {
    print(paste0("Using depth: ", depth2))
    musr2 = mus1[mus1Rnd]*depth2
  } else {
    musr2 = mus2[sample(length(mus2), numGenes, replace=TRUE)]*depth2
  }
  
  for(i in 1:numGenes) {
    lst = genGeneData(numCells, ge1[i], musr1[i])
    lst2 = genGeneData(numCells, ge2[i], musr2[i])
    geMatrix[i,1:numCells] = lst[[1]]
    geMatrix[i,((1:numCells) + numCells)] = lst2[[1]]
    histGroup1[i,] = lst[[2]]
    histGroup2[i,] = lst2[[2]]
  }


  corrMatrix = geMatrix
  
  #Predict using ZTNB
  if (!useBinDowns) {
    print("Predicting...")
    gex1 = rowSums(geMatrix[1:numGenes,1:numCells])
    gex2 = rowSums(geMatrix[1:numGenes,(1:numCells) + numCells])
    gexPred1 = rep(NA, numGenes)
    gexPred2 = rep(NA, numGenes)
    for (i in 1:numGenes) {
      gexPred1[i] = predZTNB(trimZeros(histGroup1[i,]), 10^20)
      gexPred2[i] = predZTNB(trimZeros(histGroup2[i,]), 10^20)
    }
    
    #scale to the same library size as before
    gexPred1Norm = gexPred1*sum(gex1)/sum(gexPred1)
    gexPred2Norm = gexPred2*sum(gex2)/sum(gexPred2)
    scale1 = gexPred1Norm/gex1
    scale2 = gexPred2Norm/gex2
    
    corrMatrix[1:numGenes, 1:numCells] = geMatrix[1:numGenes, 1:numCells]*scale1
    corrMatrix[1:numGenes, (1:numCells)+numCells] = geMatrix[1:numGenes, (1:numCells)+numCells]*scale2
  } else {
    print("Running binomial downsampling...")
    #use the downsampling value sent in
    depthChange = depth2/depth1
    if (depthChange > 1) { # the second group should be downsampled
      histGroup = histGroup2
      frac = 1/depthChange
    } else {
      histGroup = histGroup1
      frac = depthChange
    }
    
    
    gex1 = rowSums(geMatrix[1:numGenes,1:numCells])
    gex2 = rowSums(geMatrix[1:numGenes,(1:numCells) + numCells])
    gexPred = rep(NA, numGenes)
    for (i in 1:numGenes) {
      gexPred[i] = binomialDowns(trimZeros(histGroup[i,]), frac)
    }
    
    #scale to the same library size as before
    if (depthChange > 1) {
      gexPred2Norm = gexPred*sum(gex1)/sum(gexPred)
      scale2 = gexPred2Norm/gex2
      corrMatrix[1:numGenes, (1:numCells)+numCells] = geMatrix[1:numGenes, (1:numCells)+numCells]*scale2
    } else {
      gexPred1Norm = gexPred*sum(gex2)/sum(gexPred)
      scale1 = gexPred1Norm/gex1
      corrMatrix[1:numGenes, 1:numCells] = geMatrix[1:numGenes, 1:numCells]*scale1
    }
  }

  #do library size normalization of the two groups - this matters for example when you have different read depths
  secFactor = sum(geMatrix[1:numGenes, 1:numCells])/sum(geMatrix[1:numGenes, (1:numCells)+numCells])
  geMatrix[1:numGenes, (1:numCells)+numCells] = geMatrix[1:numGenes, (1:numCells)+numCells] * secFactor
  secFactor = sum(corrMatrix[1:numGenes, 1:numCells])/sum(corrMatrix[1:numGenes, (1:numCells)+numCells])
  corrMatrix[1:numGenes, (1:numCells)+numCells] = corrMatrix[1:numGenes, (1:numCells)+numCells] * secFactor
  
  

  #extract fold changes for false positives and false negatives
  lfcUnc = log2(rowSums(geMatrix[1:numGenes,1:numCells])/rowSums(geMatrix[1:numGenes,((1:numCells) + numCells)]))
  

  #extract fold changes for false positives and false negatives
  lfcCorr = log2(rowSums(corrMatrix[1:numGenes,1:numCells])/rowSums(corrMatrix[1:numGenes,((1:numCells) + numCells)]))

  return(list(isDE, lfcUnc, lfcCorr))
}

#Test 2 - Test LFCs for depth
#No depth difference, should give clear answers
res = evaluateCondition(standardExpr, mus1=c(1,0.25), mus2=NA, depth2 = 1, useBinDowns = TRUE, DEExprFactor = 2, numGenes=500)
hist(res[[2]]) #should produce 3 spikes, at -1, 0, and 1, where 0 is twice as high, which it does, so ok

#now, with depth difference
res = evaluateCondition(standardExpr, mus1=c(5,0.1), mus2=NA, depth2 = 0.1, useBinDowns = TRUE, DEExprFactor = 2, numGenes=500)
#calc possible outcomes
expr1 = 1-dnbinom(0, mu=5, size=1) #0.8333333
expr2 = 1-dnbinom(0, mu=0.5, size=1) #0.3333333
expr3 = 1-dnbinom(0, mu=0.1, size=1) #0.09090909
expr4 = 1-dnbinom(0, mu=0.01, size=1) #0.00990099
scaleBef = (expr1 + expr3)/2
scaleAft = (expr2 + expr4)/2
#left side of pair
normExpr5 = expr1/scaleBef #1.803279
normExpr0_1 = expr3/scaleBef #0.1967213
#right side of pair
normExpr0_5 = expr2/scaleAft #1.942308
normExpr0_01 = expr4/scaleAft #0.05769231

#it results in two cases:
LFCHighCopy = log2(normExpr5/normExpr0_5) #-0.1071494
LFCLowCopy = log2(normExpr0_1/normExpr0_01) #1.769702
#total expected: -1.11, -0.11, 0.77, 0.89, 1.77, 2.77
plot(density(res[[2]], adjust=0.1)) #looks good, spikes at the right positions

#now, with no DE (DE tested above) but different mus for left and right
res = evaluateCondition(standardExpr, mus1=c(5,0.1), mus2=c(5,0.1), useBinDowns = FALSE, DEExprFactor = 1, numGenes=500)
#expected results are 
logDEFromSat = log2(normExpr5/normExpr0_1)#3.196397
#so, we would expect 25% -3.20, 50% 0 and 25% + 3.20
hist(res[[2]])#looks reasonable





#depth
#mu
#LFC
#Gene expression, i.e. number of molecules

set.seed(1)


#############################
# Depth
#############################

#We set the depth2 parameter
depths = 2^((-5):5) #5 is about the limit before the histograms starts to get truncated because of the limited length of the histograms
depthDE = list()
depthFCUnc = list()
depthFCCorr = list()

for(i in 1:length(depths)) {
  print(paste0(i, " of ", length(depths)))
  res = evaluateCondition(standardExpr, mus1=mus, mus2=NA, depth2 = depths[i], useBinDowns = TRUE, DEExprFactor = 2, numGenes=5000)
  depthDE = c(depthDE, list(res[[1]]))
  depthFCUnc = c(depthFCUnc, list(res[[2]]))
  depthFCCorr = c(depthFCCorr, list(res[[3]]))
}

depthData = list(de = depthDE, 
                 fcUnc = depthFCUnc,
                 fcCorr = depthFCCorr)
saveRDS(depthData, file =paste0(figure_data_path, "simDepthData.RDS"))

#############################
# mu
#############################
muVals = 2^seq(-7,3,by=1)
mu1 = muVals[length(muVals)] #compare all the others against the highest value

muFCUnc = rep(NA,length(muVals))
muFCCorr = rep(NA,length(muVals))

#just generate one gene with a lot of cells, highly expressed, at two different mu:s

g1 = genGeneData(5000, 1000, mu1, totMoleculesPerCell=12000)
unc1 = sum(g1[[1]])
pred1 = predZTNB(g1[[2]], 10^20)
  
for(i in 1:length(muVals)) {
  print(paste0(i, " of ", length(muVals)))
  
  g2 = genGeneData(5000, 1000, muVals[i], totMoleculesPerCell=12000)
  unc2 = sum(g2[[1]])
  pred2 = predZTNB(g2[[2]], 10^20)
  
  muFCUnc[i] = abs(log2(unc1/unc2))
  muFCCorr[i] = abs(log2(pred1/pred2))
}

muData = list(muFCUnc = muFCUnc, 
              muFCCorr = muFCCorr)
saveRDS(muData, file =paste0(figure_data_path, "simMuData.RDS"))

#############################
# Fold change
#############################

fcs = 2^(seq(0.5,3,by=(0.5)))
fcDE = list()
fcFCUnc = list()
fcFCCorr = list()

for(i in 1:length(fcs)) {
  print(paste0(i, " of ", length(fcs)))
  res = evaluateCondition(standardExpr, mus1=mus, mus2=mus, DEExprFactor = fcs[i], numGenes=5000)
  fcDE = c(fcDE, list(res[[1]]))
  fcFCUnc = c(fcFCUnc, list(res[[2]]))
  fcFCCorr = c(fcFCCorr, list(res[[3]]))
}

fcData = list(de = fcDE,
              fcUnc = fcFCUnc,
              fcCorr = fcFCCorr)
saveRDS(fcData, file =paste0(figure_data_path, "simFcData.RDS"))

#############################
# Gene expression
#############################

gexs = 2^(2:13)
gexDE = list()
gexFCUnc = list()
gexFCCorr = list()

for(i in 1:length(gexs)) {
  print(paste0(i, " of ", length(gexs)))
  res = evaluateCondition(gexs[i], mus1=mus, mus2=mus, DEExprFactor = 2, numGenes=5000)
  gexDE = c(gexDE, list(res[[1]]))
  gexFCUnc = c(gexFCUnc, list(res[[2]]))
  gexFCCorr = c(gexFCCorr, list(res[[3]]))
}

gexData = list(de = gexDE,
               fcUnc = gexFCUnc,
               fcCorr = gexFCCorr)
saveRDS(gexData, file =paste0(figure_data_path, "simGexData.RDS"))




#estimate counts per cell
expr = rep(NA,length(mus))
for(i in 1:length(mus)) {
  expr[i] = mean(genGeneData(500, 1000, mus[i])[[1]])
}

avgGeneCounts = mean(expr)
expGeneCounts = 1000/10^6*12000

countsPerCell = 12000*avgGeneCounts/expGeneCounts #7395.409

sizeFac = 1000*length(mus)/10^6
sum(expr)


mean(genGeneData(500, 1000, 0.125)[[1]])
mean(genGeneData(500, 1000, 8)[[1]])

