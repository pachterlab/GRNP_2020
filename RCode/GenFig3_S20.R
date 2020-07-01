#before anything else, you need to setup paths:
#for example
#source("C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/paths.R")
#or, in colab:
#source("GRNP_2020/RCode/pathsGoogleColab.R")

source(paste0(sourcePath, "CCCHelpers.R"))

library(ggplot2)
library(ggpubr)
library(hexbin)
library(dplyr)


ldata = readRDS(paste0(figure_data_path, "Fig3_ldata.RDS"))
ldata2 = readRDS(paste0(figure_data_path, "Fig3_ldata2.RDS"))



#generate plot data
plotdata = tibble(gene=ldata$gene, 
                  x=ldata$x, 
                  nopred=ldata$nopred - ldata$trueval,
                  pred=ldata$pred - ldata$trueval,
                  poolpred=ldata$poolpred - ldata$trueval)

#melt
plotdata.m = reshape2::melt(plotdata, id.vars=c("gene","x"), measure.vars = c("nopred", "pred", "poolpred"))

labl = labeller(variable = 
                      c("nopred" = "No Correction",
                        "pred" = "Correction",
                        "poolpred" = "Correction using Pooling"))

dfline = data.frame(x=c(0,16), y=c(0,0))

dummyData = data.frame(x=c(0,0), y=c(1.1, -1.5)) #used in a trick to set y axis range below

fig3 = ggplot(plotdata.m) +
  stat_binhex(bins=60,na.rm = TRUE, mapping=aes(x = x, y=value, fill = log(..count..))) + # opts(aspect.ratio = 1) +
  facet_wrap(facets = ~  variable, scales = "free_x", labeller = labl, ncol=3) +
  geom_line(data=dfline, mapping = aes(x=x, y=y), color="black", size=1) + 
  geom_blank(data = dummyData, mapping = aes(x=x, y=y)) + #trick to set y axis range
  labs(y=expression(Log[2]*" fold change (CPM)"), x=expression(Log[2]*"(CPM + 1)")) +
  theme(panel.background = element_rect("white", "white", 0, 0, "white"),
        legend.position= "bottom", legend.direction = "horizontal",#, legend.title = element_blank())
        strip.text.x = element_text(size = 12, face = "bold"),
        #legend.position= "none",
        strip.background = element_blank())

fig3 # for some reason this plot sometimes fail and show an error ("hbin" ...) - Restart R and try again in that case


ggsave(
  paste0(figure_path, "Fig3.png"),
  plot = fig3, device = "png",
  width = 7, height = 4, dpi = 300)


#########################
# Fig S20 (Sampling noise)
#########################
#cpm and log transform

plotdata2 = tibble(gene=ldata2$gene, 
                  x=ldata2$x, 
                  y=ldata2$sampling - ldata$nopred)


dfline = data.frame(x=c(0,16), y=c(0,0))
dummyData = data.frame(x=c(0,0), y=c(1.1, -1.5))

figS20 = ggplot(plotdata2) +
  stat_binhex(bins=60,na.rm = TRUE, mapping=aes(x = x, y=y, fill = log(..count..))) + # opts(aspect.ratio = 1) +
  #facet_wrap(facets = ~  variable, scales = "free_x", labeller = labl, ncol=3) +
  geom_line(data=dfline, mapping = aes(x=x, y=y), color="black", size=1) + 
  geom_blank(data = dummyData, mapping = aes(x=x, y=y)) + #trick to set y axis range
  labs(y=expression(Log[2]*" fold change (CPM)"), x=expression(Log[2]*"(CPM + 1)")) +
  theme(panel.background = element_rect("white", "white", 0, 0, "white"),
        legend.position= "bottom", legend.direction = "horizontal",#, legend.title = element_blank())
        strip.text.x = element_text(size = 12, face = "bold"),
        #legend.position= "none",
        strip.background = element_blank())

figS20 # for some reason this plot sometimes fail and show an error ("hbin" ...) - Restart R and try again in that case


ggsave(
  paste0(figure_path, "FigS20.png"),
  plot = figS20, device = "png",
  width = 3, height = 4, dpi = 300)



#The data to present over the plots


print(paste0("CCC, no pred: ", getCCC(ldata$nopred, ldata$trueval))) #0.981275291888894
print(paste0("CCC, pred no pooling: ", getCCC(ldata$pred, ldata$trueval))) #0.993829998877551
print(paste0("CCC, pred with pooling: ", getCCC(ldata$poolpred, ldata$trueval))) #0.997028015743896
print(paste0("CCC, no pred, ds 100 times vs ds: ", getCCC(ldata2$nopred, ldata2$sampling))) #0.99890437562339




