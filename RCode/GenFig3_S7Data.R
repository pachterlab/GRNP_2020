sourcePath = "C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/"

source(paste0(sourcePath,"ButterflyHelpers.R"))


####################################################


dsid = "PBMC_V3_3"
otherIds = c("PBMC_V3", "PBMC_V3_2", "PBMC_NG", "PBMC_NG_2", "PBMC_V2", "EVALPBMC")
loadBugs(dsid)
bugs = get(paste0("dsBugs",dsid), envir=.GlobalEnv)

dsBug = bugs[[2]]
rm(bugs, dsBugsPBMC_V3_3  ) #save some memory


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
load(file=paste0(figure_data_path, "PBMC_V3_3_ds10_20Times.RData")) #gets the data for comparison with sampling noise

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
colnames(PBMC_V3_3_ds10_20Times) = c("gene", "sampling")



#merge all
m1 = inner_join(fromStats, pred100From10, by="gene")
m2 = inner_join(m1, predPool, by="gene")

#move sampling noise to a supporting figure
#m3 = inner_join(m2, PBMC_V3_3_ds10_20Times, by="gene")

ldata = m2

m3 = inner_join(fromStats, PBMC_V3_3_ds10_20Times, by="gene")

ldata2 = m3


#cpm and log transform
#for (i in 2:7) {
for (i in 2:6) {
  ldata[, i] = log2(ldata[, i]*10^6/sum(ldata[, i]) + 1)
}

for (i in 2:5) {
  ldata2[, i] = log2(ldata2[, i]*10^6/sum(ldata2[, i]) + 1)
}


saveRDS(ldata, paste0(figure_data_path, "Fig3_ldata.RDS"))
saveRDS(ldata2, paste0(figure_data_path, "Fig3_ldata2.RDS"))

