#before anything else, you need to setup paths:
#for example
#source("C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/paths.R")
#or, in colab:
#source("GRNP_2020/RCode/pathsGoogleColab.R")

source(paste0(sourcePath,"ButterflyHelpers.R"))


library(ggplot2)
library(ggpubr)
library(hexbin)



#####################################################
# Fig S1-S3
#####################################################

loadStats("LC")
loadStats("PBMC_NG")
loadStats("PBMC_NG_2")
loadStats("PBMC_V3")
loadStats("PBMC_V3_2")
loadStats("PBMC_V3_3")
loadStats("PBMC_V2")
loadStats("EVAL")
loadStats("EVALPBMC")
loadStats("EVALPBMC_DS")
loadStats("EVALPBMC_SW")
loadStats("MRET")
loadStats("MRET2")
loadStats("MARSSEQ")


FSCMHistPlot = function(umis, fracones, title) {
  filt = fracones[umis >= 30]
  h = hist(filt, breaks=seq(0, 1, by=0.05), plot=F)
  df = data.frame(x = h$mids, y = h$density*0.05) #convert density to fraction of the total molecules
  fig = ggplot(df,aes(x=x,y=y)) +
    geom_bar(stat="identity", orientation = "x") + #orientation x is a workaround for a bug in ggplot
    labs(y="Frac. of mol.", x="FSCM", title=title) +
    theme(panel.background = element_rect("white", "white", 0, 0, "white")) +
    theme(axis.title = element_text(size = 10))
  print(fig)
  v = var(filt)
  print(paste0("Variance: ", v))
  return(list(f=fig, v= v))
}

#Fig S1

figS1_1 = FSCMHistPlot(statsPBMC_V3$UMIs_PBMC_V3_d_100, statsPBMC_V3$FracOnes_PBMC_V3_d_100, "PBMC_V3")
figS1_2 = FSCMHistPlot(statsPBMC_V3_2$UMIs_PBMC_V3_2_d_100, statsPBMC_V3_2$FracOnes_PBMC_V3_2_d_100, "PBMC_V3_2")
figS1_3 = FSCMHistPlot(statsPBMC_V3_3$UMIs_PBMC_V3_3_d_100, statsPBMC_V3_3$FracOnes_PBMC_V3_3_d_100, "PBMC_V3_3")
figS1_4 = FSCMHistPlot(statsPBMC_NG$UMIs_PBMC_NG_d_100, statsPBMC_NG$FracOnes_PBMC_NG_d_100, "PBMC_NG")
figS1_5 = FSCMHistPlot(statsPBMC_NG_2$UMIs_PBMC_NG_2_d_100, statsPBMC_NG_2$FracOnes_PBMC_NG_2_d_100, "PBMC_NG_2")
figS1_6 = FSCMHistPlot(statsPBMC_V2$UMIs_PBMC_V2_d_100, statsPBMC_V2$FracOnes_PBMC_V2_d_100, "PBMC_V2")
figS1_7 = FSCMHistPlot(statsEVAL$UMIs_EVAL_d_100, statsEVAL$FracOnes_EVAL_d_100, "EVAL")
figS1_8 = FSCMHistPlot(statsEVALPBMC$UMIs_EVALPBMC_d_100, statsEVALPBMC$FracOnes_EVALPBMC_d_100, "EVALPBMC")
figS1_9 = FSCMHistPlot(statsLC$UMIs_LC_d_100, statsLC$FracOnes_LC_d_100, "LC")
figS1_10 = FSCMHistPlot(statsMRET2$UMIs_MRET2_d_100, statsMRET2$FracOnes_MRET2_d_100, "MRET2")

figS1 = ggarrange(figS1_1$f, figS1_2$f, figS1_3$f, figS1_4$f, figS1_5$f, figS1_6$f, figS1_7$f,
                       figS1_8$f, figS1_9$f, figS1_10$f, nrow=4, ncol=3,
                       labels=c("A","B","C","D","E","F","G","H","I","J"))
figS1
ggsave(
  paste0(figure_path, "FigS1.png"),
  plot = figS1, device = "png",
  width = 6, height = 7, dpi = 300)

tenxvariances = c(figS1_1$v,figS1_2$v,figS1_3$v,figS1_4$v,figS1_5$v,figS1_6$v,figS1_7$v,figS1_8$v,figS1_9$v,figS1_10$v)
mean(tenxvariances)#0.03326073
sd(tenxvariances)#0.01171921
#hist(tenxvariances, 10)


figS2_1 = FSCMHistPlot(statsEVALPBMC_DS$UMIs_EVALPBMC_DS_d_100, statsEVALPBMC_DS$FracOnes_EVALPBMC_DS_d_100, "EVALPBMC_DS")
figS2_2 = FSCMHistPlot(statsMRET$UMIs_MRET_d_100, statsMRET$FracOnes_MRET_d_100, "MRET")

figS2 = ggarrange(figS2_1$f, figS2_2$f, nrow=1, ncol=2,
                  labels=c("A","B"))
figS2
ggsave(
  paste0(figure_path, "FigS2.png"),
  plot = figS2, device = "png",
  width = 4, height = 2, dpi = 300)

dsvariances = c(figS2_1$v,figS2_2$v)
mean(dsvariances)#0.01790591
sd(dsvariances)#0.001654036

figS3_1 = FSCMHistPlot(statsEVALPBMC_SW$UMIs_EVALPBMC_SW_d_100, statsEVALPBMC_SW$FracOnes_EVALPBMC_SW_d_100, "EVALPBMC_SW")
figS3_2 = FSCMHistPlot(statsMARSSEQ$UMIs_MARSSEQ_d_100, statsMARSSEQ$FracOnes_MARSSEQ_d_100, "MARSSEQ")

figS3 = ggarrange(figS3_1$f, figS3_2$f, nrow=1, ncol=2,
                  labels=c("A","B"))

#figS3 = ggarrange(figS3_1, nrow=1, ncol=1,
#                       labels=c("A"))
figS3
ggsave(
  paste0(figure_path, "FigS3.png"),
  plot = figS3, device = "png",
  width = 4, height = 2, dpi = 300)

swvariances = c(figS3_1$v)
mean(swvariances)#0.02304077
MARSSEQvariances = c(figS3_2$v)
mean(MARSSEQvariances)#0.01653287



