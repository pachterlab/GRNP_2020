#
# Generates the data for the plots in Fig. 1 and 3
#

#before anything else, you need to setup paths:
#for example
#source("C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/paths.R")
#or, in colab:
#source("GRNP_2020/RCode/pathsGoogleColab.R")

source(paste0(sourcePath,"ButterflyHelpers.R"))
source(paste0(sourcePath,"PreseqHelpers.R"))

loadStats("EVAL")
#so, use the histograms from downsampled data at 0.25, which somewhat matches the A figure
loadBug("EVAL", 0.25)

#Fig 3A - histograms per gene 


#collapsedNonFilt = bugEVAL %>% group_by(gene) %>% do(countslist=c(.$count))
collapsedNonFilt = bug_EVAL_25 %>% group_by(gene) %>% do(countslist=c(.$count))

h1 = hist(collapsedNonFilt$countslist[collapsedNonFilt$gene == "Vmn1r13"][[1]], breaks=seq(0.5, 100.5, by=1), plot=F)
h2 = hist(collapsedNonFilt$countslist[collapsedNonFilt$gene == "Ubb"][[1]], breaks=seq(0.5, 100.5, by=1), plot=F)

saveRDS(h1, paste0(figure_data_path, "Fig3_h1.RDS"))
saveRDS(h2, paste0(figure_data_path, "Fig3_h2.RDS"))


#now, fig 3C

#create prediction data

xes = c(1,2,4,5,8,12,16,20)
predVals = 20/xes
downSamp = c(0.05, 0.1, 0.2, 0.25, 0.4, 0.6, 0.8, 1)

#build it backwards
cpms = tibble(gene=statsEVAL$gene, n=statsEVAL$CPM_EVAL_d_100)
for (i in (length(xes)-1):1) {
  loadBug("EVAL", downSamp[i])
  pred = upSampleAndGetMeanExprPreSeqZTNB(getBug("EVAL", downSamp[i]), t=predVals[[i]])
  rmBug("EVAL", downSamp[i])
  cpm = pred
  cpm[[2]] = cpm[[2]]*10^6/sum(cpm[[2]])
  cpms = inner_join(cpm, cpms, by="gene")
}

r1 = data.frame(cpms)[cpms$gene == "Vmn1r13",]
r2 = data.frame(cpms)[cpms$gene == "Ubb",]

saveRDS(r1, paste0(figure_data_path, "Fig3C_r1.RDS"))
saveRDS(r2, paste0(figure_data_path, "Fig3C_r2.RDS"))
