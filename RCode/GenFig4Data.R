sourcePath = "C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/"

source(paste0(sourcePath,"ButterflyHelpers.R"))

library(ggplot2)
library("ggpubr")

#prediction is slow, so make it possible to do this once and save it
genFig4Data = function(dsid1, dsid2, predVals) {
  bug1 = get(paste0("bug",dsid1), envir=.GlobalEnv)
  bug2 = get(paste0("bug",dsid2), envir=.GlobalEnv)
  #predict and add results to merged
  pred1 = upSampleAndGetMeanExprPreSeqZTNB(bug1, t=predVals)
  pred2 = upSampleAndGetMeanExprPreSeqZTNB(bug2, t=predVals)
  #pred1 = upSampleAndGetMeanExprPreSeq(bug1, t=predVals, mt=2)
  #pred2 = upSampleAndGetMeanExprPreSeq(bug2, t=predVals, mt=2)
  return(list(pred1,pred2))
}


#Two plots; EVALPBMC_DS vs EVALPBMC and EVALPBMC_SW vs EVALPBMC

predVals_1 = c(1, 1.5, 2, 3, 4, 6, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096)

loadBug("EVALPBMC_DS")
loadBug("EVALPBMC")

d1 = genFig4Data("EVALPBMC_DS", "EVALPBMC", predVals_1)


saveRDS(d1, paste0(figure_data_path, "Fig4_d1.RDS"))


