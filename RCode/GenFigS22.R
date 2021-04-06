#before anything else, you need to setup paths:
#for example
#source("C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/paths.R")
#or, in colab:
#source("GRNP_2020/RCode/pathsGoogleColab.R")

source(paste0(sourcePath,"ButterflyHelpers.R"))
source(paste0(sourcePath,"CCCHelpers.R"))

library(ggplot2)
library("ggpubr")


###########################################
#Load the data
###########################################

depthData = readRDS(file=paste0(figure_data_path, "simDepthData.RDS"))
muData = readRDS(file =paste0(figure_data_path, "simMuData.RDS"))
fcData = readRDS(file =paste0(figure_data_path, "simFcData.RDS"))
gexData = readRDS(file =paste0(figure_data_path, "simGexData.RDS"))

numGenes = 500
logDepths = (-5):5
logMus = seq(-7,3,by=1)
lfcs = seq(0.5,3,by=0.5)
logGexs = 2:13


#install.packages("pROC")
library(pROC)

####################
# FC
####################

ns = length(fcData$de)

x=lfcs[1:ns]
y1=rep(NA,ns)
y2=rep(NA,ns)

for (i in 1:ns) {
  rocInfoUnc = roc(fcData$de[[i]], abs(fcData$fcUnc[[i]]), plot=FALSE )
  rocInfoCorr = roc(fcData$de[[i]], abs(fcData$fcCorr[[i]]), plot=FALSE )
  y1[i] = rocInfoUnc$auc
  y2[i] = rocInfoCorr$auc
}

group = factor(c(rep(0,ns),rep(1,ns)), c(0,1), c("Uncorrected", "Corrected"))

dsPlot = data.frame(x=c(x,x), y=c(y1,y2), group=group)

pA = ggplot(dsPlot, aes(x=x, y=y, group=group, color=group)) + 
  geom_line() + #geom_point(size=1.5) +
  labs(y="AUC", x="DE log fold change") +
  ggtitle("DE LFC") +
  theme(panel.background = element_rect("white", "white", 0, 0, "white"), legend.position= "bottom", 
        legend.direction = "horizontal", legend.box = "horizontal", legend.title = element_blank()) +
  guides(fill=guide_legend(nrow=1, ncol=2, byrow=F))

print(pA)

####################
# Depth
####################

ns = length(depthData$de)

x=logDepths
y1=rep(NA,ns)
y2=rep(NA,ns)

for (i in 1:ns) {
    rocInfoUnc = roc(depthData$de[[i]], abs(depthData$fcUnc[[i]]), plot=FALSE )
    rocInfoCorr = roc(depthData$de[[i]], abs(depthData$fcCorr[[i]]), plot=FALSE )
    y1[i] = rocInfoUnc$auc
    y2[i] = rocInfoCorr$auc
}

group = factor(c(rep(0,ns),rep(1,ns)), c(0,1), c("Uncorrected", "Corrected"))

dsPlot = data.frame(x=c(x,x), y=c(y1,y2), group=group)

pB = ggplot(dsPlot, aes(x=x, y=y, group=group, color=group)) + 
  geom_line() + #geom_point(size=1.5) +
  labs(y="AUC", x=expression(Log[2]*"(rel. read depth)")) +
  ggtitle("Read depth") +
  theme(panel.background = element_rect("white", "white", 0, 0, "white"), legend.position= "bottom", 
        legend.direction = "horizontal", legend.box = "horizontal", legend.title = element_blank()) +
  guides(fill=guide_legend(nrow=1, ncol=2, byrow=F))

print(pB)




####################
# Mu
####################

x=logMus

y1=muData$muFCUnc
y2=muData$muFCCorr

group = factor(c(rep(0,ns),rep(1,ns)), c(0,1), c("Uncorrected", "Corrected"))

dsPlot = data.frame(x=c(x,x), y=c(y1,y2), group=group)

pC = ggplot(dsPlot, aes(x=x, y=y, group=group, color=group)) + 
  geom_line() + #geom_point(size=1.5) +
  labs(y="Log fold change", x="NB mean") +
  ggtitle("NB size LFC") +
  theme(panel.background = element_rect("white", "white", 0, 0, "white"), legend.position= "bottom", 
        legend.direction = "horizontal", legend.box = "horizontal", legend.title = element_blank()) +
  guides(fill=guide_legend(nrow=1, ncol=2, byrow=F))

print(pC)


####################
# Gene expr
####################

ns = length(gexData$de)

x=logGexs
y1=rep(NA,ns)
y2=rep(NA,ns)

for (i in 1:ns) {
  rocInfoUnc = roc(gexData$de[[i]], abs(gexData$fcUnc[[i]]), plot=FALSE )
  rocInfoCorr = roc(gexData$de[[i]], abs(gexData$fcCorr[[i]]), plot=FALSE )
  y1[i] = rocInfoUnc$auc
  y2[i] = rocInfoCorr$auc
}

group = factor(c(rep(0,ns),rep(1,ns)), c(0,1), c("Uncorrected", "Corrected"))

dsPlot = data.frame(x=c(x,x), y=c(y1,y2), group=group)

pD = ggplot(dsPlot, aes(x=x, y=y, group=group, color=group)) + 
  geom_line() + #geom_point(size=1.5) +
  labs(y="AUC", x=expression(Log[2]*"(CPM)")) +
  ggtitle("Gene expression") +
  theme(panel.background = element_rect("white", "white", 0, 0, "white"), legend.position= "bottom", 
        legend.direction = "horizontal", legend.box = "horizontal", legend.title = element_blank()) +
  guides(fill=guide_legend(nrow=1, ncol=2, byrow=F))

print(pD)


###########################
# Assemble all plots
###########################

figSim = ggarrange(pC, pA, pB, pD, nrow=2, ncol=2, labels=c("A","B","C","D"))

ggsave(
  paste0(figure_path, "FigS22.png"),
  plot = figSim, device = "png",
  width = 6, height = 6, dpi = 300)

