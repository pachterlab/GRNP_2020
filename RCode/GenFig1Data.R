#
# Generates the data for the plots in Fig. 1
#
sourcePath = "C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/"

source(paste0(sourcePath,"ButterflyHelpers.R"))

loadBugs("EVAL")
loadStats("EVAL")
loadBug("EVAL")

#Fig 1B, I - histograms per gene 

#so, create the histograms from downsampled data at 0.25, which somewhat matches the A figure
load(file=paste0(figure_data_path, "EVAL/bugEVAL025.RData"))

#collapsedNonFilt = bugEVAL %>% group_by(gene) %>% do(countslist=c(.$count))
collapsedNonFilt = bugEVAL025 %>% group_by(gene) %>% do(countslist=c(.$count))

h1 = hist(collapsedNonFilt$countslist[collapsedNonFilt$gene == "Vmn1r13"][[1]], breaks=seq(0.5, 100.5, by=1), plot=F)
h2 = hist(collapsedNonFilt$countslist[collapsedNonFilt$gene == "Ubb"][[1]], breaks=seq(0.5, 100.5, by=1), plot=F)

saveRDS(h1, paste0(figure_data_path, "Fig1_h1.RDS"))
saveRDS(h2, paste0(figure_data_path, "Fig1_h2.RDS"))


#now, fig 1 B II and III


#fig1B_III
#create prediction data

xes = c(1,2,4,8,12,16,20)
predVals = 20/xes

#build it backwards
cpms = tibble(gene=statsEVAL$gene, n=statsEVAL$CPM_EVAL_d_100)
for (i in (length(xes)-1):1) {
  pred = upSampleAndGetMeanExprPreSeqZTNB(dsBugsEVAL[[i]], t=predVals[[i]])
  cpm = pred
  cpm[[2]] = cpm[[2]]*10^6/sum(cpm[[2]])
  cpms = inner_join(cpm, cpms, by="gene")
}

r1 = data.frame(cpms)[cpms$gene == "Vmn1r13",]
r2 = data.frame(cpms)[cpms$gene == "Ubb",]

saveRDS(r1, paste0(figure_data_path, "Fig1_r1_III.RDS"))
saveRDS(r2, paste0(figure_data_path, "Fig1_r2_III.RDS"))


pred025 = upSampleAndGetMeanExprPreSeqZTNB(bugEVAL025, t=4)
pred025$CPM = pred025$p4 * 10^6/sum(pred025$p4)
r1_025 = pred025$CPM[pred025$gene == "Vmn1r13"]
r2_025 = pred025$CPM[pred025$gene == "Ubb"]

saveRDS(r1_025, paste0(figure_data_path, "Fig1_r1_025_III.RDS"))
saveRDS(r2_025, paste0(figure_data_path, "Fig1_r2_025_III.RDS"))


