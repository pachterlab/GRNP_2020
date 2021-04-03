#before anything else, you need to setup paths:
#for example
#source("C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/paths.R")
#or, in colab:
#source("GRNP_2020/RCode/pathsGoogleColab.R")

source(paste0(sourcePath,"ButterflyHelpers.R"))
source(paste0(sourcePath,"preseqHelpers.R"))

library(tictoc)

GenAlgEvaluationData <- function(dsid) {
  loadStats(dsid)
  stats = getStats(dsid)
  
  
  ##################################
  #first, prediction of all UMIs
  ##################################
  
  #Get UMI counts at all different ds stages
  ###########################################
  
  dss = c("5","10","20","40","60","80","100")
  
  #generate the columns to extract from stats
  colList = dss #allocation of right size
  for (i in 1:length(dss)) {
    colList[i] = paste0("UMIs_", dsid, "_d_", dss[i])
  }
  extrFromStats = stats[,which(colnames(stats) %in% colList)]
  dsCounts = as.numeric(colSums(extrFromStats))

  
  #get histogram from 0.05:
  loadBug(dsid, 0.05)
  h = totalCPUHistogram(getBug(dsid, 0.05))
  rmBug(dsid, 0.05)

  x = c(0.05,0.1,0.2,0.4,0.6,0.8,1)
  t = c(1,2,4,8,12,16,20)


  #predict with Good-Toulmin
  predGT005 = rep(0,length(t))
  predGT005[1] = dsCounts[1];
  for (i in 2:length(t)) {
    predGT005[i] = goodToulmin(h,t[i])
  }
  
  #predict with Preseq DS, mt=20
  predPSDS005_20 = rep(0,length(t))
  predPSDS005_20[1] = dsCounts[1];
  for (i in 2:length(t)) {
    predPSDS005_20[i] = predPreSeqDS(h,t[i],20)
  }
  
  #predict with Preseq DS, mt=2
  predPSDS005_2 = rep(0,length(t))
  predPSDS005_2[1] = dsCounts[1];
  for (i in 2:length(t)) {
    predPSDS005_2[i] = predPreSeqDS(h,t[i],2)
  }
  
  #predict with Preseq ZTNB
  predPSZTNB005 = rep(0,length(t))
  predPSZTNB005[1] = dsCounts[1];
  for (i in 2:length(t)) {
    predPSZTNB005[i] = predPreSeqZTNB(h,t[i]) #just ignore the warnings, deprecated...
  }
  
  ##################################
  #now, prediction per gene for different methods
  ##################################
  
  loadBug(dsid, 0.1)
  ds10Bug = getBug(dsid, 0.1)
  rmBug(dsid, 0.1)

  collapsed10 = aggregate(count~gene, ds10Bug, FUN=c) #if you get an error here, you probably defined a variable called "c"...
  totUMIs10 = sapply(collapsed10$count, FUN=length)
  rm(ds10Bug)
  
  loadBug(dsid, 1)
  ds100Bug = getBug(dsid, 1)
  rmBug(dsid, 1)
  
  collapsedFull = aggregate(count~gene, ds100Bug, FUN=c) #if you get an error here, you probably defined a variable called "c"...
  rm(ds100Bug)
  merged2 = inner_join(collapsed10, collapsedFull, by="gene")
  rm(collapsed10, collapsedFull)
  colnames(merged2) = c("gene", "DS10", "Full")
  
  totUMIsFull = sapply(merged2$Full, FUN=length)
  #sort the genes on number of UMIs
  srt = sort(totUMIs10, index.return=T)
  umis = srt$x
  merged2srt = merged2[srt$ix,]
  
  #now predict using both ztnb and ds:
  numgenes = dim(merged2srt)[1]
  predds_20 = rep(0,numgenes)
  predds_2 = rep(0,numgenes)
  predztnb = rep(0,numgenes)
  predbp = rep(0,numgenes)
  predscaled = rep(0,numgenes)
  fracOnes = rep(0,numgenes)
  fullUMIs = rep(0,numgenes)
  
  globScale = sum(totUMIsFull)/sum(totUMIs10)

  #All quick ones.
  print(paste0("DS etc.:", numgenes))
  for (i in 1:numgenes) {
    if (i %% 1000 == 0) {
      print(i)
    }
    h = hist(merged2srt$DS10[[i]], breaks=seq(0.5, max(merged2srt$DS10[[i]])+0.5, by=1), plot = F)
    predds_20[[i]] = predPreSeqDS(h, 10, 20)
    predds_2[[i]] = predPreSeqDS(h, 10, 2)
    predscaled[[i]] = length(merged2srt$DS10[[i]])*globScale #this resembles CPM
    
    fullUMIs[[i]] = length(merged2srt$Full[[i]])
    
    fracOnes[[i]] = h$density[1]
  }
  
  #ztnb
  print(paste0("ZTNB:", numgenes))
  for (i in 1:numgenes) {
    if (i %% 1000 == 0) {
      print(i)
    }
    h = hist(merged2srt$DS10[[i]], breaks=seq(0.5, max(merged2srt$DS10[[i]])+0.5, by=1), plot = F)
    predztnb[[i]] = predPreSeqZTNB(h, 10)
  }
  
  #best practice
  print(paste0("Best practice:", numgenes))
  for (i in 1:numgenes) {
    if (i %% 500 == 0) {
      print(i)
    }
    h = hist(merged2srt$DS10[[i]], breaks=seq(0.5, max(merged2srt$DS10[[i]])+0.5, by=1), plot = F)
    predbp[[i]] = predPreSeq(h, 10, mt=2)
  }
  
  
  toSave = list(x, dsCounts, predGT005, predPSDS005_2, predPSDS005_20, predPSZTNB005, 
                predds_20, predds_2, predztnb, predbp, predscaled, fracOnes, fullUMIs, umis,
                merged2srt)
  
  
  filename = paste0(figure_data_path, dsid, "/PredEvalData.RDS")
  saveRDS(toSave, filename)
  
  
}

GenAlgEvaluationData("EVAL")
GenAlgEvaluationData("EVALPBMC")
GenAlgEvaluationData("EVALPBMC_DS")
GenAlgEvaluationData("EVALPBMC_SW")
GenAlgEvaluationData("PBMC_V3")
GenAlgEvaluationData("PBMC_V3_2")
GenAlgEvaluationData("PBMC_V3_3")
GenAlgEvaluationData("PBMC_NG")
GenAlgEvaluationData("PBMC_NG_2")
GenAlgEvaluationData("PBMC_V2")
GenAlgEvaluationData("LC")
GenAlgEvaluationData("MRET")
GenAlgEvaluationData("MRET2")
GenAlgEvaluationData("MARSSEQ")
