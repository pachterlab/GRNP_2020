#before anything else, you need to setup paths:
#for example
#source("C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/paths.R")
#or, in colab:
#source("GRNP_2020/RCode/pathsGoogleColab.R")

source(paste0(sourcePath,"ButterflyHelpers.R"))

library(ggplot2)
library(ggpubr)
library(hexbin)
library(tidyverse)


#####################################################
# Fig S4
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


AddToHexbinData = function(dat, umis, fracOnes, dataset) {
  logUmis = log2(umis)
  d = tibble(x = logUmis, y = fracOnes, ds = rep(dataset, length(fracOnes)))

  if (is.null(dat)) {
    dat = d
  } else {
    dat = bind_rows(dat,d);
  }

  return(dat)
}


dat = NULL

dat = AddToHexbinData(dat, statsPBMC_V3$UMIs_PBMC_V3_d_100, statsPBMC_V3$FracOnes_PBMC_V3_d_100, "PBMC_V3")
dat = AddToHexbinData(dat, statsPBMC_V3_2$UMIs_PBMC_V3_2_d_100, statsPBMC_V3_2$FracOnes_PBMC_V3_2_d_100, "PBMC_V3_2")
#dat = AddToHexbinData(dat, statsPBMC_V3_3$UMIs_PBMC_V3_3_d_100, statsPBMC_V3_3$FracOnes_PBMC_V3_3_d_100, "PBMC_V3_3") #this is the same as Fig 2A, don't duplicate it!
dat = AddToHexbinData(dat, statsPBMC_NG$UMIs_PBMC_NG_d_100, statsPBMC_NG$FracOnes_PBMC_NG_d_100, "PBMC_NG")
dat = AddToHexbinData(dat, statsPBMC_NG_2$UMIs_PBMC_NG_2_d_100, statsPBMC_NG_2$FracOnes_PBMC_NG_2_d_100, "PBMC_NG_2")
dat = AddToHexbinData(dat, statsPBMC_V2$UMIs_PBMC_V2_d_100, statsPBMC_V2$FracOnes_PBMC_V2_d_100, "PBMC_V2")
dat = AddToHexbinData(dat, statsEVAL$UMIs_EVAL_d_100, statsEVAL$FracOnes_EVAL_d_100, "EVAL")
dat = AddToHexbinData(dat, statsEVALPBMC$UMIs_EVALPBMC_d_100, statsEVALPBMC$FracOnes_EVALPBMC_d_100, "EVALPBMC")
dat = AddToHexbinData(dat, statsLC$UMIs_LC_d_100, statsLC$FracOnes_LC_d_100, "LC")
dat = AddToHexbinData(dat, statsMRET2$UMIs_MRET2_d_100, statsMRET2$FracOnes_MRET2_d_100, "MRET2")
dat = AddToHexbinData(dat, statsEVALPBMC_DS$UMIs_EVALPBMC_DS_d_100, statsEVALPBMC_DS$FracOnes_EVALPBMC_DS_d_100, "EVALPBMC_DS")
dat = AddToHexbinData(dat, statsMRET$UMIs_MRET_d_100, statsMRET$FracOnes_MRET_d_100, "MRET")
dat = AddToHexbinData(dat, statsEVALPBMC_SW$UMIs_EVALPBMC_SW_d_100, statsEVALPBMC_SW$FracOnes_EVALPBMC_SW_d_100, "EVALPBMC_SW")
dat = AddToHexbinData(dat, statsMARSSEQ$UMIs_MARSSEQ_d_100, statsMARSSEQ$FracOnes_MARSSEQ_d_100, "MARSSEQ")

dat$ds = factor(dat$ds, levels = c("PBMC_V3","PBMC_V3_2","PBMC_V3_3","PBMC_NG","PBMC_NG_2","PBMC_V2","EVAL","EVALPBMC","LC","MRET2","EVALPBMC_DS","MRET","EVALPBMC_SW","MARSSEQ"))

#create figure:
figS4 = ggplot(dat) +
  stat_binhex(bins=60,na.rm = TRUE, mapping=aes(x = x, y=y, fill = log(..count..))) + # opts(aspect.ratio = 1) +
  facet_wrap(facets = ~ds, scales = "free_x", ncol=3) +
  labs(x=expression(Log[2]*"(UMI counts)"), y="FSCM") +
  theme(panel.background = element_rect("white", "white", 0, 0, "white"),
        legend.position= "bottom", legend.direction = "horizontal",#, legend.title = element_blank())
        strip.text.x = element_text(size = 12, face = "bold"),
        #legend.position= "none",
        strip.background = element_blank())

figS4 # for some reason this plot sometimes fail and show an error ("hbin" ...) - Restart R and try again in that case

ggsave(
  paste0(figure_path, "FigS4.png"),
  plot = figS4, device = "png",
  width = 7, height = 11, dpi = 300)


#############################
# Fig S5 - Hexbin plot
#############################

AddToHexbinData2 = function(dat, ds1, ds2) {
  stats1 = get(paste0("stats",ds1), envir=.GlobalEnv)
  indUMIs1 = which(colnames(stats1) == paste0("UMIs_",ds1,"_d_100"))
  indFO1 = which(colnames(stats1) == paste0("FracOnes_",ds1,"_d_100"))
  stats2 = get(paste0("stats",ds2), envir=.GlobalEnv)
  indUMIs2 = which(colnames(stats2) == paste0("UMIs_",ds2,"_d_100"))
  indFO2 = which(colnames(stats2) == paste0("FracOnes_",ds2,"_d_100"))
  
  stats1Filt = stats1[stats1[[indUMIs1]] >= 200, ]
  stats2Filt = stats2[stats2[[indUMIs2]] >= 200, ]
  
  merged = inner_join(stats2Filt[, c(1,indFO2)], stats1Filt[, c(1,indFO1)], by="gene")
  colnames(merged) = c("gene", "x", "y")
  
  d = merged %>% add_column(ds=paste0(ds1, " vs ", ds2))
  

  if (is.null(dat)) {
    dat = d
  } else {
    dat = bind_rows(dat,d);
  }
  
  return(dat)
}

#ds1 = "PBMC_V3_3"
#ds2 = "PBMC_V3_2"

dat2 = NULL

#dat2 = AddToHexbinData2(dat2, "PBMC_V3_3", "PBMC_V3_2") #the same as Fig 2B, don't duplicate it!
dat2 = AddToHexbinData2(dat2, "PBMC_V3_3", "PBMC_V2")
dat2 = AddToHexbinData2(dat2, "PBMC_V2", "EVALPBMC")
dat2 = AddToHexbinData2(dat2, "PBMC_V2", "LC")
dat2 = AddToHexbinData2(dat2, "EVALPBMC", "EVALPBMC_DS")
dat2 = AddToHexbinData2(dat2, "EVALPBMC", "EVALPBMC_SW")
dat2 = AddToHexbinData2(dat2, "EVALPBMC_DS", "EVALPBMC_SW")
dat2 = AddToHexbinData2(dat2, "MRET2", "MRET")
dat2 = AddToHexbinData2(dat2, "EVAL", "MARSSEQ")

#dat2 = AddToHexbinData2(dat2, "PBMC_V3_3", "PBMC_V2")

#specify the order of the plots
dat2$ds = factor(dat2$ds, levels = c("PBMC_V3_3 vs PBMC_V3_2", "PBMC_V3_3 vs PBMC_V2", "PBMC_V2 vs EVALPBMC",
                                     "PBMC_V2 vs LC", "EVALPBMC vs EVALPBMC_DS", "EVALPBMC vs EVALPBMC_SW", 
                                     "EVALPBMC_DS vs EVALPBMC_SW", "MRET2 vs MRET", "EVAL vs MARSSEQ"))

dfline = data.frame(x=c(0,1), y=c(0,1))

figS5 = ggplot(dat2) +
  stat_binhex(bins=60,na.rm = TRUE, mapping=aes(x = x, y=y, fill = log(..count..))) + # opts(aspect.ratio = 1) +
  geom_line(data=dfline, mapping=aes(x = x, y=y), color="black", size=1.5) + 
  facet_wrap(facets = ~ds, scales = "free_x", ncol=2) +
  labs(x="FSCM 2", y="FSCM 1") +
  theme(panel.background = element_rect("white", "white", 0, 0, "white"),
        legend.position= "bottom", legend.direction = "horizontal",#, legend.title = element_blank())
        strip.text.x = element_text(size = 10, face = "bold"),
        #legend.position= "none",
        strip.background = element_blank())

figS5 # for some reason this plot sometimes fail and show an error ("hbin" ...) - Restart R and try again in that case

ggsave(
  paste0(figure_path, "FigS5.png"),
  plot = figS5, device = "png",
  width = 6, height = 10, dpi = 300)


#############################
# Fig 2 - One of each type above
#############################

dat3 = NULL
dat3 = AddToHexbinData(dat3, statsPBMC_V3_3$UMIs_PBMC_V3_3_d_100, statsPBMC_V3_3$FracOnes_PBMC_V3_3_d_100, "PBMC_V3_3")

fig2A = ggplot(dat3) +
  stat_binhex(bins=60,na.rm = TRUE, mapping=aes(x = x, y=y, fill = log(..count..))) + # opts(aspect.ratio = 1) +
  #facet_wrap(facets = ~ds, scales = "free_x", ncol=3) +
  ggtitle("FSCM vs Gene Expression") + 
  labs(x=expression(Log[2]*"(UMI counts)"), y="FSCM") +
  theme(panel.background = element_rect("white", "white", 0, 0, "white"),
        legend.position= "bottom", legend.direction = "horizontal",#, legend.title = element_blank())
        strip.text.x = element_text(size = 12, face = "bold"),
        #legend.position= "none",
        strip.background = element_blank())

fig2A

dat4 = NULL
dat4 = AddToHexbinData2(dat4, "PBMC_V3_3", "PBMC_V3_2")
dfline = data.frame(x=c(0,1), y=c(0,1))

fig2B = ggplot(dat4) +
  stat_binhex(bins=60,na.rm = TRUE, mapping=aes(x = x, y=y, fill = log(..count..))) + # opts(aspect.ratio = 1) +
  geom_line(data=dfline, mapping=aes(x = x, y=y), color="black", size=1.5) + 
  ggtitle("FSCM Across Datasets") + 
  #facet_wrap(facets = ~ds, scales = "free_x", ncol=2) +
  labs(x="FSCM, PBMC_V3_2", y="FSCM, PBMC_V3_3") +
  theme(panel.background = element_rect("white", "white", 0, 0, "white"),
        legend.position= "bottom", legend.direction = "horizontal",#, legend.title = element_blank())
        strip.text.x = element_text(size = 10, face = "bold"),
        #legend.position= "none",
        strip.background = element_blank())

fig2B # for some reason this plot sometimes fail and show an error ("hbin" ...) - Restart R and try again in that case

fig2 = ggarrange(fig2A, fig2B, nrow=1, ncol=2,labels=c("A","B"))

fig2 # for some reason this plot sometimes fail and show an error ("hbin" ...) - Restart R and try again in that case

ggsave(
  paste0(figure_path, "Fig2.png"),
  plot = fig2, device = "png",
  width = 6, height = 4, dpi = 300)



