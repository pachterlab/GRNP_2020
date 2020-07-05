#before anything else, you need to setup paths:
#for example
#source("C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/paths.R")
#or, in colab:
#source("GRNP_2020/RCode/pathsGoogleColab.R")

source(paste0(sourcePath,"ButterflyHelpers.R"))
source(paste0(sourcePath,"CCCHelpers.R"))

library(ggplot2)
library("ggpubr")


###########################################
#Assemble the data
###########################################


dsid1 = "EVALPBMC_DS"
dsid2 = "EVALPBMC"

predVals = c(1, 1.5, 2, 3, 4, 6, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096)
title = "Drop-Seq vs Chromium v2"

loadStats("EVALPBMC_DS")
loadStats("EVALPBMC")

#get precalculated data (generated using GenFig4)
precalcData = readRDS(paste0(figure_data_path, "Fig4_d1.RDS"))

#l1 = genFig4Plot("EVALPBMC_DS", "EVALPBMC", predVals_1, 100, "Drop-Seq vs Chromium v2", d1)

cpmfiltLimit = 100

#genFig4Plot <- function(dsid1, dsid2, predVals, cpmfiltLimit, title, precalcData) {
stats1 = get(paste0("stats",dsid1), envir=.GlobalEnv)
stats2 = get(paste0("stats",dsid2), envir=.GlobalEnv)

pred1 = precalcData[[1]];
pred2 = precalcData[[2]];

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

merged1 = inner_join(ds1, pred1, by="gene")
merged2 = inner_join(ds2, pred2, by="gene")
#make sure that prediction has not changed p1 (i.e. predicting to 1 times the counts)
merged1$p1 = merged1$CPM_EVALPBMC_DS_d_100
merged2$p1 = merged2$CPM_EVALPBMC_d_100

mergedTmp = inner_join(merged1, merged2, by="gene")
for (i in 2:(dim(mergedTmp)[2])) {
  mergedTmp[,i] = mergedTmp[,i]*10^6/sum(mergedTmp[,i])
}

#filter out the lowly expressed genes - prediction is too uncertain for those genes to show any interesting results
mergedTot = mergedTmp[(mergedTmp$p1.x > cpmfiltLimit) & (mergedTmp$p1.y > cpmfiltLimit),]

#mergeTot first has all columns of dataset 1 in order, then the columns of dataset 2

print(paste0("Number of genes: ", dim(mergedTot)[1]))

#now CPM again to remove the uncertainty in the lowly expressed genes and log transform
for (i in 2:(dim(mergedTot)[2])) {
  mergedTot[,i] = log2(mergedTot[,i]*10^6/sum(mergedTot[,i]) + 1)
}

xes = c(log2(dsx), log2(predVals))


###########################################
#Fig A
###########################################


#The plot has 2 curves, both just 10x prediction and prediction of both datasets
dsLen = length(dsGrade)
predLen = length(predVals)
numPoints = dsLen + predLen
corsBothPred = rep(0,numPoints)
for (i in 1:numPoints) {
  corsBothPred[i] = getCCC(mergedTot[[i+1]], mergedTot[[i+numPoints+1]])
}

#now create the plot with 2 curves, both double prediction and just 10x prediction
pred = factor(c(rep(0,dsLen),rep(1, length(predVals))), c(0,1), c("Downsampled", "Predicted"))
x = xes
y = corsBothPred

df = data.frame(x=x, y=y, pred = pred)

fig4A = ggplot(df, aes(x=x, y=y, group=pred, linetype=pred)) + 
  geom_line() + #geom_point(size=1.5) +
  labs(y="CCC", x=expression(Log[2]*"(Prediction range)")) +
  ggtitle("Concordance Between Datasets") +
  theme(panel.background = element_rect("white", "white", 0, 0, "white"), legend.position= "bottom", 
        legend.direction = "vertical", legend.box = "horizontal", legend.title = element_blank()) +
  guides(fill=guide_legend(nrow=2, ncol=1, byrow=F))

print(fig4A)

#save them individually since we want different heights, so the x axis labels of both plots are at the same y coord:
ggsave(
  paste0(figure_path, "Fig4A.png"),
  plot = fig4A, device = "png",
  width = 3.65, height = 4.3, dpi = 300)


###########################################
#Fig 4B
###########################################

g1 = "MTRNR2L1"
g2 = "LSP1"
g3 = "EMP3"

gene1 = as.numeric(mergedTot[mergedTot$gene == g1, -1])
gene2 = as.numeric(mergedTot[mergedTot$gene == g2, -1])
gene3 = as.numeric(mergedTot[mergedTot$gene == g3, -1])

genes = c(gene1,gene2,gene3)
geneNames = c(g1,g2,g3)

#create data frame for plot:
y = genes
x = rep(xes,2*length(geneNames))
pred = rep(factor(c(rep(0,dsLen),rep(1, length(predVals))), c(0,1), c("Downs.", "Pred.")),2*length(geneNames))
gene = rep(0:(length(geneNames)-1), each=length(xes)*2)
genef = factor(gene, 0:(length(geneNames)-1), geneNames)
technology = factor(rep(c(rep(0,numPoints), rep(1,numPoints)),length(geneNames)), c(0,1), c("Drop-Seq", "10X"))

df = data.frame(x=x, y=y, dspred = pred, Gene = genef, Technology = technology, stringsAsFactors = F)

fig4B = ggplot(df, aes(x=x, y=y, colour = Gene, linetype=dspred, shape=Technology, group=interaction(Gene,dspred, Technology))) + 
  geom_line() + geom_point(size=1.2) +
  ggtitle("Prediction per Gene") +
  xlab(expression(Log[2]*"(sampling factor)")) + ylab(expression(Log[2]*"(CPM + 1)")) +
  theme(panel.background = element_rect("white", "white", 0, 0, "white"),
        legend.position= "bottom", legend.direction = "vertical", legend.box = "horizontal", legend.title = element_blank(),
        strip.background = element_blank()) +
  guides(fill=guide_legend(nrow=2, ncol=1, byrow=F))
fig4B


ggsave(
  paste0(figure_path, "Fig4B.png"),
  plot = fig4B, device = "png",
  width = 3.65, height = 4.545, dpi = 300)




#some test code to find interesting genes
##########################

#get stats
#mergedStats = inner_join(stats1, stats2, by="gene")
#filtStats = mergedStats$CPM_EVALPBMC_DS_d_100 > 100 & mergedStats$CPM_EVALPBMC_d_100 > 100
#mergedStatsFilt1 = mergedStats[filtStats,]
#dim(mergedStatsFilt1)

#selStats = abs(mergedStatsFilt1$FracOnes_EVALPBMC_DS_d_100 - mergedStatsFilt1$FracOnes_EVALPBMC_d_100) > 0.4
#mergedStatsFilt1$gene[selStats]

#diffs = abs(mergedStatsFilt1$FracOnes_EVALPBMC_DS_d_100 - mergedStatsFilt1$FracOnes_EVALPBMC_d_100)

#selStats = diffs > 0.2 & diffs < 0.22 & mergedTot$CPM_EVALPBMC_DS_d_100 < 120
#selStats = diffs < 0.05 & mergedStatsFilt1$CPM_EVALPBMC_DS_d_100 > 500
#improvement = abs(mergedTot$CPM_EVALPBMC_DS_d_100 - mergedTot$CPM_EVALPBMC_d_100) - abs(mergedTot$p1024.x - mergedTot$p1024.y)
#selStats = improvement > 2
#mergedStatsFilt1$gene[selStats]




