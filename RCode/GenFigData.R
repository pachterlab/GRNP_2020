sourcePath = "C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/"

source(paste0(sourcePath,"ButterflyHelpers.R"))
source(paste0(sourcePath,"GenBugSummary.R"))


#Generate BUG data from the datasets
#########################################

createStandardBugsData(paste0(dataPath,"EVAL/"), "EVAL", c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1))
#createStandardBugsData(paste0(dataPath,"EVAL_NUC/"), "EVAL_NUC", c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1))
createStandardBugsData(paste0(dataPath,"EVALPBMC/"), "EVALPBMC", c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1))
createStandardBugsData(paste0(dataPath,"EVALPBMCDropseq/"), "EVALPBMC_DS", c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1), 200)
#createStandardBugsData(paste0(dataPath,"EVALPBMC_CELSeq2/"), "EVALPBMC_CS", c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1), 200)
createStandardBugsData(paste0(dataPath,"EVALPBMC_SeqWell/"), "EVALPBMC_SW", c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1), 200)
createStandardBugsData(paste0(dataPath,"lc/"), "LC", c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1), 1000)#large dataset, need to discard at a different level here
createStandardBugsData(paste0(dataPath,"pbmc5kNextGEM/"), "PBMC_NG", c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1))
createStandardBugsData(paste0(dataPath,"pbmcv3/"), "PBMC_V3", c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1))
createStandardBugsData(paste0(dataPath,"pbmcv3_2/"), "PBMC_V3_2", c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1))
createStandardBugsData(paste0(dataPath,"pbmcv3_3/"), "PBMC_V3_3", c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1))
createStandardBugsData(paste0(dataPath,"pbmcv2/"), "PBMC_V2", c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1))
createStandardBugsData(paste0(dataPath,"PBMC_NG2/"), "PBMC_NG_2", c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1))
#createStandardBugsData(paste0(dataPath,"FISH_CMP/"), "FISH_CMP", c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1))
createStandardBugsData(paste0(dataPath,"MRET/"), "MRET", c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1))
createStandardBugsData(paste0(dataPath,"MRET2/"), "MRET2", c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1))


#figure out suitable human genes for many and few single-copy molecules

loadBug("EVALPBMC")

#find the 100 most highly expressed genes
collapsed = aggregate(count~gene, bugEVALPBMC, FUN=c) #if you get an error here, you probably defined a variable called "c"...

counts = bugEVALPBMC %>% group_by(gene) %>% tally()
srt = sort(counts$n, index.return=T, decreasing=T)
nn = 200
genes = counts$gene[srt$ix[1:nn]]
colsel = collapsed[collapsed$gene %in% genes,]
fracOnes = rep(0,nn)
for (i in 1:nn) {
  fracOnes[i] = sum(colsel$count[[i]] == 1)/length(colsel$count[[i]])
}
plot(fracOnes)
sel = fracOnes > 0.8
colsel$gene[sel]#FGF23 looks good for high fraction of ones


hist(bugEVALPBMC[bugEVALPBMC$gene == "FGF23",]$count, breaks=seq(0.5, max(bug$count)+0.5, by=1))#looks good

sel = fracOnes < 0.01
colsel$gene[sel]#RPS10 looks good for high fraction of ones
hist(bugEVALPBMC[bugEVALPBMC$gene == "RPS10",]$count, breaks=seq(0.5, max(bug$count)+0.5, by=1))#looks good

#create statistics
genBugSummary("EVAL", "Vmn1r13", "Ubb", 10)
genBugSummary("EVALPBMC", "FGF23", "RPS10", 10)
genBugSummary("EVALPBMC_DS", "FGF23", "RPS10", 10)
genBugSummary("EVALPBMC_CS", "FGF23", "RPS10", 10)
genBugSummary("EVALPBMC_SW", "FGF23", "RPS10", 10)
genBugSummary("LC", "FGF23", "RPS10", 10)
genBugSummary("PBMC_NG", "FGF23", "RPS10", 12)
genBugSummary("PBMC_V3", "FGF23", "RPS10", 12)
genBugSummary("PBMC_V3_2", "FGF23", "RPS10", 12)
genBugSummary("PBMC_V3_3", "FGF23", "RPS10", 12)
genBugSummary("PBMC_V2", "FGF23", "RPS10", 10)
genBugSummary("PBMC_NG_2", "FGF23", "RPS10", 12)
genBugSummary("FISH_CMP", "FGF23", "RPS10", 12)
genBugSummary("MRET", "Vmn1r13", "Ubb", 10)
genBugSummary("MRET2", "Vmn1r13", "Ubb", 10)


loadBug("PBMC_V3")
loadBug("PBMC_V3_2")
loadBug("PBMC_V3_3")
loadBug("PBMC_NG")
loadBug("PBMC_NG_2")
loadBug("PBMC_V2")
loadBug("EVALPBMC")
#generate the pooled histograms
generatePooledHistogramDS("PBMC_V3")
generatePooledHistogramDS("PBMC_V3_2")
generatePooledHistogramDS("PBMC_V3_3")
generatePooledHistogramDS("PBMC_NG")
generatePooledHistogramDS("PBMC_NG_2")
generatePooledHistogramDS("PBMC_V2")
generatePooledHistogramDS("EVALPBMC")

#test that it worked
loadPooledHistogramDS("PBMC_V3")
pooledHistDSPBMC_V3[1:10,1:10] #looks reasonable

loadPooledHistogramDS("PBMC_V3_2")
loadPooledHistogramDS("PBMC_V3_3")
loadPooledHistogramDS("PBMC_NG")
loadPooledHistogramDS("PBMC_NG_2")
loadPooledHistogramDS("PBMC_V2")
loadPooledHistogramDS("EVALPBMC")



generatePooledHistogram("PBMC_V3")
generatePooledHistogram("PBMC_V3_2")
generatePooledHistogram("PBMC_V3_3")
generatePooledHistogram("PBMC_NG")
generatePooledHistogram("PBMC_NG_2")
generatePooledHistogram("PBMC_V2")
generatePooledHistogram("EVALPBMC")

#test that it worked
loadPooledHistogram("PBMC_V3")
pooledHistPBMC_V3[1:10,1:10] #looks reasonable



#generate downsampling 20 times of PBMCV3_3 - to enable comparison with sampling noise
loadBug("PBMC_V3_3")
PBMC_V3_3_ds10_20Times = downSampleBUGNTimes(bugPBMC_V3_3, 0.1, 20)
save(PBMC_V3_3_ds10_20Times, file=paste0(figure_data_path, "PBMC_V3_3_ds10_20Times.RData"))

#generate downsampling to 0.3 for EVAL - for fig 1
loadBug("EVAL")
bugEVAL025 = downSampleBUG(bugEVAL, 0.25)
save(bugEVAL025, file=paste0(figure_data_path, "EVAL/bugEVAL025.RData"))

#Generate data for multimapped figure:
######################################

bugFileEVAL = read.table(paste0(dataPath,"EVAL/bus_output/bug.txt"), stringsAsFactors = F)


#split up the "gene1,gene2" in the fourth column into lists
#geneLists = lapply(bugFileEVAL[,4], function(s) strsplit(s, ",")[[1]])
#gllenghts = sapply(geneLists, length)

colnames(bugFileEVAL) = c("barcode", "umi", "gene", "count")

library(stringr)
gllenghts = sapply(bugFileEVAL[,3],function(s) str_count(s, fixed(",")) + 1)
names(gllenghts) = NULL

#Get multimapped reads

bugEVAL_MultiMapped = bugFileEVAL[gllenghts > 1,]

save(bugEVAL_MultiMapped, file=paste0(figure_data_path, "EVAL/bugMulti.RData"))


