#before anything else, you need to setup paths:
#for example
#source("C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/paths.R")
#or, in colab:
#source("GRNP_2020/RCode/pathsGoogleColab.R")

source(paste0(sourcePath,"ButterflyHelpers.R"))
source(paste0(sourcePath,"preseqHelpers.R"))
source(paste0(sourcePath,"BinomialDownsampling.R"))


####################################################


dsid = "PBMC_V3_3"
otherIds = c("PBMC_V3", "PBMC_V3_2", "PBMC_NG", "PBMC_NG_2", "PBMC_V2", "EVALPBMC")
loadBug(dsid, 0.1)

dsBug = getBug(dsid, 0.1)


loadPooledHistogramDS("PBMC_V3_3")
loadPooledHistogramDS("PBMC_V3_2")
loadPooledHistogramDS("PBMC_V3")
loadPooledHistogramDS("PBMC_NG")
loadPooledHistogramDS("PBMC_NG_2")
loadPooledHistogramDS("PBMC_V2")
loadPooledHistogramDS("EVALPBMC")



#Collect the data

poolHistList = poolHistograms(dsid, dsBug, otherIds)

loadStats(dsid)

#create data for supplementary plot
loadBug(dsid)
bug = getBug(dsid)
binDs = binomialDownsampling(bug, 0.1)




#no prediction
fromStats = tibble(gene = statsPBMC_V3_3$gene, 
                   trueval = statsPBMC_V3_3$CPM_PBMC_V3_3_d_100,
                   x = statsPBMC_V3_3$CPM_PBMC_V3_3_d_10, #so, we use the 
                   nopred = statsPBMC_V3_3$CPM_PBMC_V3_3_d_10)

#prediction
#pred100From10 = upSampleAndGetMeanExprPreSeq(dsBug, t=10, mt=2)
pred100From10 = upSampleAndGetMeanExprPreSeqZTNB(dsBug, t=10)

colnames(pred100From10) = c("gene", "pred")

#prediction with pooling
predPool = poolPrediction(dsBug, 10, poolHistList, 500000)
colnames(predPool) = c("gene", "poolpred")



#sampling noise
colnames(binDs) = c("gene", "sampling")



#merge all
m1 = inner_join(fromStats, pred100From10, by="gene")
m2 = inner_join(m1, predPool, by="gene")

#move sampling noise to a supporting figure

ldata = m2

m3 = inner_join(fromStats, binDs, by="gene")

ldata2 = m3


#cpm and log transform
#for (i in 2:7) {
for (i in 2:6) {
  ldata[, i] = log2(ldata[, i]*10^6/sum(ldata[, i]) + 1)
}

for (i in 2:5) {
  ldata2[, i] = log2(ldata2[, i]*10^6/sum(ldata2[, i]) + 1)
}


saveRDS(ldata, paste0(figure_data_path, "Fig4AC_ldata.RDS"))
saveRDS(ldata2, paste0(figure_data_path, "Fig4AC_ldata2.RDS"))

