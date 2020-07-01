#before anything else, you need to setup paths:
#for example
#source("C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/paths.R")
#or, in colab:
#source("GRNP_2020/RCode/pathsGoogleColab.R")

source(paste0(sourcePath,"ButterflyHelpers.R"))
source(paste0(sourcePath,"CCCHelpers.R"))

library(ggplot2)
library("ggpubr")

#assumes the stats are loaded
#expects the first in predVals to be 1
genFig4Plot <- function(dsid1, dsid2, predVals, cpmfiltLimit, title, precalcData) {
  stats1 = get(paste0("stats",dsid1), envir=.GlobalEnv)
  stats2 = get(paste0("stats",dsid2), envir=.GlobalEnv)

  pred1 = precalcData[[1]];
  pred2 = precalcData[[2]];
  
  
  merged = inner_join(pred1, pred2, by="gene")
  #CPM for filtering below
  for (i in 2:(2*length(predVals)+1)) {
    merged[,i] = merged[,i]*10^6/sum(merged[,i])
  }
  
  #extract down-sampled expression from stats:
  dsx = c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1)
  dsGrade = c(5,10,20,40,60,80,100)
  ids1 = NULL
  ids2 = NULL
  for (gr in dsGrade) {
    ids1 = c(ids1, paste0("CPM_", dsid1, "_d_", gr))
    ids2 = c(ids2, paste0("CPM_", dsid2, "_d_", gr))
  }
  ds1 = stats1[,c(1,which(colnames(stats1) %in% ids1))]
  ds2 = stats2[,c(1,which(colnames(stats2) %in% ids2))]
  mergedDs = inner_join(ds1, ds2, by="gene")
 
  #filter out the lowly expressed genes - prediction is too uncertain for those genes to show any interesting results
  filt = merged[(merged$p1.x > cpmfiltLimit) & (merged$p1.y > cpmfiltLimit),]

  #Now merge predicted and downsampled
  mergedTot = inner_join(filt, mergedDs, by="gene")
  
  print(paste0("Number of genes: ", dim(mergedTot)[1]))
  
  #now CPM again to remove the uncertainty in the lowly expressed genes
  for (i in 2:(dim(mergedTot)[2])) {
    mergedTot[,i] = log2(mergedTot[,i]*10^6/sum(mergedTot[,i]) + 1)
  }
  
  #now create the plot
  cors = rep(0,length(predVals))
  for (i in 1:(length(predVals))) {
    cors[i] = getCCC(mergedTot[,i+1, drop=T], mergedTot[,i+length(predVals)+1, drop=T])
  }
  
  
  dsLen = length(dsGrade)
  cors2 = rep(0,dsLen)
  dsOffs =  2*length(predVals)+1
  for (i in 1:dsLen) {
    cors2[i] = getCCC(mergedTot[,dsOffs + i, drop=T], mergedTot[,dsOffs + i + dsLen, drop=T])
  }
  
  pred = factor(c(rep(0,dsLen),rep(1, length(predVals))), c(0,1), c("Downsampled", "Predicted"))
  
  df = data.frame(x=c(log2(dsx), log2(predVals)), y=c(cors2,cors), pred = pred)

  
  fig = ggplot(df, aes(x=x, y=y, linetype= pred)) + 
    geom_line() + geom_point() +
    labs(y=bquote(R^2~"vs Identity"), x=expression(Log[2]*"(Prediction range)"), title = title) +
    theme(panel.background = element_rect("white", "white", 0, 0, "white"), legend.position= "none")
  
  figLegend = ggplot(df, aes(x=x, y=y, linetype= pred)) + 
    geom_line() + geom_point() +
    labs(y="CCC", x=expression(Log[2]*"(Prediction range)")) +
    theme(panel.background = element_rect("white", "white", 0, 0, "white"), legend.position= "bottom", 
          legend.direction = "horizontal", legend.title = element_blank())
  
  print(fig)
  print(figLegend)
  
  return(list(fig, figLegend))
}


predVals_1 = c(1, 1.5, 2, 3, 4, 6, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096)
title = "10x Chromium v2 vs Drop-Seq"

loadStats("EVALPBMC_DS")
loadStats("EVALPBMC")

#get precalculated data (generated using GenFig4)
d1 = readRDS(paste0(figure_data_path, "Fig4_d1.RDS"))

l1 = genFig4Plot("EVALPBMC_DS", "EVALPBMC", predVals_1, 100, "Drop-Seq vs Chromium v2", d1)

fig4 = l1[[2]]

ggsave(
  paste0(figure_path, "Fig4.png"),
  plot = fig4, device = "png",
  width = 3, height = 3, dpi = 300)




