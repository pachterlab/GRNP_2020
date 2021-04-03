#before anything else, you need to setup paths:
#for example
#source("C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/paths.R")
#or, in colab:
#source("GRNP_2020/RCode/pathsGoogleColab.R")

library(ggplot2)
library("ggpubr")


###########################################
#Load the data
###########################################

gctlcu = readRDS(file=paste0(figure_data_path, "gc.RDS"))

###########################################
#Create the plots
###########################################
dsPlot = tibble(x=gctlcu$gc, y=gctlcu$FSCM)
loess_fit = loess(y~x, dsPlot, span = 0.3)
dsLoess = data.frame(x = dsPlot$x, y = predict(loess_fit))

pA = ggplot(dsPlot, aes(x=x, y=y)) +
  stat_binhex(bins=60,na.rm = TRUE, mapping=aes(x = x, y=y, fill = log(..count..))) + # opts(aspect.ratio = 1) +
  geom_line(data=dsLoess, color="#FF8800", size=1.4) + 
  labs(y="FSCM", x="GC content", title="FSCM vs GC content") +
  theme(panel.background = element_rect("white", "white", 0, 0, "white"),
        legend.position= "bottom", legend.direction = "horizontal",#, legend.title = element_blank())
        strip.text.x = element_text(size = 12, face = "bold"),
        #legend.position= "none",
        strip.background = element_blank())
pA

dsPlot = tibble(x=gctlcu$gc, y=gctlcu$CU)
loess_fit = loess(y~x, dsPlot, span = 0.3)
dsLoess = data.frame(x = dsPlot$x, y = predict(loess_fit))

pB = ggplot(dsPlot, aes(x=x, y=y)) +
  stat_binhex(bins=60,na.rm = TRUE, mapping=aes(x = x, y=y, fill = log(..count..))) + # opts(aspect.ratio = 1) +
  geom_line(data=dsLoess, color="#FF8800", size=1.4) + 
  labs(y="Copies per UMI", x="GC content", title="CU vs GC content") +
  theme(panel.background = element_rect("white", "white", 0, 0, "white"),
        legend.position= "bottom", legend.direction = "horizontal",#, legend.title = element_blank())
        strip.text.x = element_text(size = 12, face = "bold"),
        #legend.position= "none",
        strip.background = element_blank())
pB

dsPlot = tibble(x=log2(gctlcu$txlen), y=gctlcu$FSCM)
dim(dsPlot)#11826
#filter extreme values
dsPlot = dsPlot[(dsPlot$x < 15) & (dsPlot$x > 7),]
dim(dsPlot)#11823, so, filtered 3 genes

loess_fit = loess(y~x, dsPlot, span = 0.3)
dsLoess = data.frame(x = dsPlot$x, y = predict(loess_fit))

pC = ggplot(dsPlot, aes(x=x, y=y)) +
  stat_binhex(bins=60,na.rm = TRUE, mapping=aes(x = x, y=y, fill = log(..count..))) + # opts(aspect.ratio = 1) +
  geom_line(data=dsLoess, color="#FF8800", size=1.4) + 
  labs(y="FSCM", x=expression(Log[2]*"(Avg. transcript length)"), title="FSCM vs Transcript Length") +
  theme(panel.background = element_rect("white", "white", 0, 0, "white"),
        legend.position= "bottom", legend.direction = "horizontal",#, legend.title = element_blank())
        strip.text.x = element_text(size = 12, face = "bold"),
        #legend.position= "none",
        strip.background = element_blank())
pC

dsPlot = tibble(x=log2(gctlcu$txlen), y=gctlcu$CU)
dim(dsPlot)#11826
#filter extreme values
dsPlot = dsPlot[(dsPlot$x < 15) & (dsPlot$x > 7),]
dim(dsPlot)#11823, so, filtered 3 genes

loess_fit = loess(y~x, dsPlot, span = 0.3)
dsLoess = data.frame(x = dsPlot$x, y = predict(loess_fit))

pD = ggplot(dsPlot, aes(x=x, y=y)) +
  stat_binhex(bins=60,na.rm = TRUE, mapping=aes(x = x, y=y, fill = log(..count..))) + # opts(aspect.ratio = 1) +
  geom_line(data=dsLoess, color="#FF8800", size=1.4) + 
  labs(y="Copies per UMI", x=expression(Log[2]*"(Avg. transcript length)"), title="CU vs Transcript Length") +
  theme(panel.background = element_rect("white", "white", 0, 0, "white"),
        legend.position= "bottom", legend.direction = "horizontal",#, legend.title = element_blank())
        strip.text.x = element_text(size = 12, face = "bold"),
        #legend.position= "none",
        strip.background = element_blank())
pD



###########################
# Assemble all plots
###########################

figGC = ggarrange(pA, pB, pC, pD, nrow=2, ncol=2, labels=c("A","B","C","D"))

ggsave(
  paste0(figure_path, "FigGC.png"),
  plot = figGC, device = "png",
  width = 6, height = 7.4, dpi = 300)

