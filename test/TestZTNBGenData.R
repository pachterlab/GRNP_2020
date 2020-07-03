#
# TCR0004 - performance increase vs accuracy loss for ZTNB
#

sourcePath = "C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/"
source(paste0(sourcePath,"ButterflyHelpers.R"))

dsid = "PBMC_V3_3"
loadBugs(dsid)
bugs = get(paste0("dsBugs",dsid), envir=.GlobalEnv)

dsBug = bugs[[2]]
rm(bugs, dsBugsPBMC_V3_3  ) #save some memory


#Collect the data

loadStats(dsid)

#no prediction
fromStats = tibble(gene = statsPBMC_V3_3$gene, 
                   trueval = statsPBMC_V3_3$CPM_PBMC_V3_3_d_100,
                   x = statsPBMC_V3_3$CPM_PBMC_V3_3_d_10, #so, we use the 
                   nopred = statsPBMC_V3_3$CPM_PBMC_V3_3_d_10)

#prediction
pred100From10Slow = upSampleAndGetMeanExprPreSeqZTNB(dsBug, t=10, incTol = 1e-10, iterIncTol = 70000)
pred100From10Fast = upSampleAndGetMeanExprPreSeqZTNB(dsBug, t=10)

colnames(pred100From10Slow) = c("gene", "predSlow")
colnames(pred100From10Fast) = c("gene", "predFast")


#merge all
m1 = inner_join(fromStats, pred100From10Slow, by="gene")
m2 = inner_join(m1, pred100From10Fast, by="gene")

testZTNBData = m2


#cpm and log transform
for (i in 2:6) {
  testZTNBData[, i] = log2(testZTNBData[, i]*10^6/sum(testZTNBData[, i]) + 1)
}

saveRDS(testZTNBData, paste0(sourcePath, "test/tmp/testZTNB.RDS"))

