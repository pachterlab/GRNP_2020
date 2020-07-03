#
# TCR0004 - performance increase vs accuracy loss for ZTNB
#

sourcePath = "C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/"
source(paste0(sourcePath,"ButterflyHelpers.R"))

library(ggplot2)
library(ggpubr)
library(hexbin)

testZTNBData = readRDS(paste0(sourcePath, "test/tmp/testZTNB.RDS"))

plotdata = tibble(gene=testZTNBData$gene, 
                  x=testZTNBData$x, 
                  nopred=testZTNBData$nopred - testZTNBData$trueval,
                  predSlow=testZTNBData$predSlow - testZTNBData$trueval,
                  predFast=testZTNBData$predFast - testZTNBData$trueval)

#melt
plotdata.m = reshape2::melt(plotdata, id.vars=c("gene","x"), measure.vars = c("nopred", "predSlow", "predFast"))

labl = labeller(variable = 
                  c("nopred" = "No Correction",
                    "predSlow" = "Correction Slow",
                    "predFast" = "Correction Fast"))

dfline = data.frame(x=c(0,16), y=c(0,0))


figTestZTNB = ggplot(plotdata.m) +
  stat_binhex(bins=60,na.rm = TRUE, mapping=aes(x = x, y=value, fill = log(..count..))) + # opts(aspect.ratio = 1) +
  facet_wrap(facets = ~  variable, scales = "free_x", labeller = labl, ncol=3) +
  geom_line(data=dfline, mapping = aes(x=x, y=y), color="black", size=1) + 
  labs(y=expression(Log[2]*" fold change (CPM)"), x=expression(Log[2]*"(CPM + 1)")) +
  theme(panel.background = element_rect("white", "white", 0, 0, "white"),
        legend.position= "bottom", legend.direction = "horizontal",#, legend.title = element_blank())
        strip.text.x = element_text(size = 12, face = "bold"),
        #legend.position= "none",
        strip.background = element_blank())

figTestZTNB # for some reason this plot sometimes fail and show an error ("hbin" ...) - Restart R and try again in that case


ggsave(
  paste0(figure_path, "FigTestZTNB.png"),
  plot = figTestZTNB, device = "png",
  width = 7, height = 4, dpi = 300)
