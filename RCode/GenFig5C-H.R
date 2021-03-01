#before anything else, you need to setup paths:
#for example
#source("C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/paths.R")
#or, 
#source("C:/Work/R/ButterflyPublicationRepo/GRNP_2020/RCode/pathsPublRepo.R")
#or, in colab:
#source("GRNP_2020/RCode/pathsGoogleColab.R")

source(paste0(sourcePath,"ButterflyHelpers.R"))
source(paste0(sourcePath,"CCCHelpers.R"))
source(paste0(sourcePath,"preseqHelpers.R"))
sourcePath2 = "C:/Work/R/ButterflyPublicationRepo/GRNP_2020/RCode/"
source(paste0(sourcePath,"ggplotHelpers.R"))


library(ggplot2)
library("ggpubr")
library(tidyverse)
library(BUSpaRse)
library(textTinyR)
library(qdapTools)
library("matrixStats")
library(Matrix.utils)

colors = gg_color_hue(2)


########################
# Load data
########################


#Seurat cannot handle the tx data
#So, cluster on gene data instead, it doesn't really matter
gn10x = read_count_output(dir=paste0(dataPath, "EVALPBMC/bus_output/genecounts"), "output", tcc=FALSE)
countsPerCell = sparse_Sums(gn10x)
gn10x = gn10x[,countsPerCell > 200]
gn10x = gn10x[sparse_Sums(gn10x, rowSums = TRUE) != 0,]
#convert the genes to gene symbols
tr2g = read.table(paste0(dataPath, "EVALPBMC/bus_output/transcripts_to_genes.txt"), stringsAsFactors = F)
lookupTable = tr2g[,2:3]
lookupTable= unique(lookupTable)

outGenes = lookup(row.names(gn10x), lookupTable)

row.names(gn10x) = outGenes
gn10x = aggregate.Matrix(gn10x, outGenes, fun='sum')



set.seed(1)

library(Seurat)

#run through Seurat
d = CreateSeuratObject(counts = gn10x, project = "10x", min.cells = 0, min.features = 0)

#create variable for mitochondrial count fraction
d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "^MT-")
VlnPlot(d, features = c("percent.mt"))

cellFilter = d$percent.mt < 10

#easisest to recreate the Seurat object

gn10x = gn10x[,cellFilter]
gn10xRawCounts = gn10xRawCounts[,cellFilter]
d = CreateSeuratObject(counts = gn10x, project = "10x", min.cells = 0, min.features = 0)


#library size normalization + change to log scale
d <- NormalizeData(d, normalization.method = "LogNormalize", scale.factor = 10000)

#Finds the most variable genes
d <- FindVariableFeatures(d, selection.method = "vst", nfeatures = 2000)


#Make mean and variance the same for all genes:
d <- ScaleData(d, features = rownames(d))

#Principal component analysis
d <- RunPCA(d, features = VariableFeatures(object = d))

#do clustering
d <- FindNeighbors(d, dims = 1:10)
d <- FindClusters(d, resolution = 0.5)

#Generate UMAP map
d <- RunUMAP(d, dims = 1:10)

clust = d$seurat_clusters

d2 = subset(d, seurat_clusters == 0 | seurat_clusters == 1 | seurat_clusters == 2 | seurat_clusters == 6  | seurat_clusters == 7 )
DimPlot(d2, reduction = "umap")
d2 = subset(d, seurat_clusters == 3 | seurat_clusters == 5 )
DimPlot(d2, reduction = "umap")
d2 = subset(d, seurat_clusters == 4)
DimPlot(d2, reduction = "umap")

sum(clust == 7)#check number of cells in cluster 7: 397

#d = RenameIdents(object = d, `0` = "T1", `1` = "T2", `2` = "T3", `3` = "M1", `4` = "B", `5` = "T4", `6` = "M2", `7` = "U1", `8` = "U2", `9` = "U3", `10` = "U4",  `11` = "U5")
d = RenameIdents(object = d, `0` = "T1", `1` = "T2", `2` = "T3", `3` = "M1", `4` = "B", `5` = "M2", `6` = "T4", `7` = "T5", `8` = "U1", `9` = "U2", `10` = "U3",  `11` = "U4")


#show clustering
pSupA = DimPlot(d, reduction = "umap")
pSupA

#FeaturePlot(d, features = c("ENST00000416293.7","ENST00000261733.6","ENST00000548536.1","ENST00000549106.1","ENST00000551450.1"))

ggsave(
  paste0(figure_path, "FigS_Seurat_5C-H.png"),
  plot = pSupA, device = "png",
  width = 4, height = 4, dpi = 300)





#calculate counts per UMI



######################
# Compare CU across clusters
######################

cuPerCluster = matrix(, nrow = nrow(gn10x), ncol = 12) #number of clusters is 12
rownames(cuPerCluster) = rownames(gn10x)
UMIsPerCluster = matrix(, nrow = nrow(gn10x), ncol = 12) #number of clusters is 12
rownames(UMIsPerCluster) = rownames(gn10x)
cellsPerCluster = rep(NA,12)
clusterNames = c("T1","T2","T3","M1","B","M2","T4","T5","U1","U2","U3","U4") 
colnames(cuPerCluster) = clusterNames
colnames(UMIsPerCluster) = clusterNames

#load the bug
loadBug("EVALPBMC")
bug = getBug("EVALPBMC")

for (i in 1:12) {
  sel = d$seurat_clusters == i-1
  cellsPerCluster[i] = sum(sel)
  UMIsPerCluster[,i] = sparse_Sums(gn10x[,sel], rowSums = TRUE)

  #to figure out cu per cluster, we need to use the bug
  barcodes = colnames(gn10x)[sel]
  subBug = bug[bug$barcode %in% barcodes,]
  
  res = subBug %>% group_by(gene) %>% summarize(cu = mean(count))
  colnames(res) = c("x","y")
  cuPerCluster[,i] = lookup(rownames(cuPerCluster), res)
}
#test ()
#sel = d$seurat_clusters == 0
#barcodes = colnames(gn10x)[sel]#814 long, reasonable
#subBug = bug[(bug$barcode %in% barcodes) & (bug$gene == "AACS"),]
#cuPerCluster[rownames(cuPerCluster) == "AACS" ,1] #2.393939
#mean(subBug$count) #2.393939 - looks good!



cuPerClusterFilt = cuPerCluster[,cellsPerCluster > 300]
UMIsPerClusterFilt = UMIsPerCluster[,cellsPerCluster > 300]
cuPerClusterFilt[UMIsPerClusterFilt < 20] = NA

#only measure the variance for somewhat highly expressed genes
cuPerClusterFilt2 = cuPerClusterFilt[(rowMeans(UMIsPerClusterFilt) > 60) & (rowSums(!is.na(cuPerClusterFilt)) > 2),]
dim(cuPerClusterFilt2)
cuPerClusterFilt2
variances = rowVars(cuPerClusterFilt2, na.rm=TRUE)

srt = sort(variances, index.return=TRUE, decreasing = TRUE)
sel = seq_len(length(srt$x)) %in% srt$ix[1:50] #get the 50 genes with most variance

#sel = variances > 10
sum(sel)

TPMs = UMIsPerCluster
for (i in 1:12) {
  TPMs[,i] = TPMs[,i]*10^6/sum(TPMs[,i])
}
colSums(TPMs)#ok


cuForVariableGenes = cuPerClusterFilt2[sel,]
UMIsForVariableGenes = UMIsPerCluster[rownames(UMIsPerCluster) %in% rownames(cuForVariableGenes),]
TPMsForVariableGenes = TPMs[rownames(TPMs) %in% rownames(cuForVariableGenes),cellsPerCluster > 300]
logtrans = log2(TPMsForVariableGenes + 1)

max(variances, na.rm=TRUE)

geneSel = rownames(cuForVariableGenes) == "ALDH2" #ALDH2
cu = cuForVariableGenes[geneSel,]
logExpr = logtrans[geneSel,]
texts = names(cu)

dsPlot = tibble(x=cu, y=logExpr, col=c(colors[1],colors[1],colors[1],colors[2],"#000000",colors[2],colors[1],colors[1]))
dsText = tibble(x=cu + 0.2, y=logExpr)
fit = lm(logExpr~cu)
summary(fit) #p = 0.0040, F-Test
#fix the texts a bit so they are visible and don't overlap
dsText$x[6] = dsText$x[6] - 1.5
dsText$y[8] = dsText$y[8] + 0.05
dsText$y[3] = dsText$y[3] - 0.05

pE = ggplot2::ggplot(dsPlot,ggplot2::aes(x=x,y=y)) +
  ggplot2::geom_point(ggplot2::aes(x=x,y=y), color=dsPlot$col, shape=1, size=1) +
  geom_abline(slope = fit$coefficients[2], intercept = fit$coefficients[1], colour="#008800", size=1.3) +
  geom_text(data=dsText, label=texts, size=4, hjust = 0, color=dsPlot$col, parse=FALSE) +
  ggplot2::labs(y=expression(Log[2]*"(CPM + 1)"), x="Copies per UMI", title="ALDH2 acr. clusters, uncorr.") +
  theme(panel.background = element_rect("white", "white", 0, 0, "white"),
        legend.position= "bottom", legend.direction = "horizontal",#, legend.title = element_blank())
        strip.text.x = element_text(size = 12, face = "bold"),
        #legend.position= "none",
        plot.title = element_text(face = "bold"),
        strip.background = element_blank())
pE




###########################
#look at the histograms for T cells and Monocytes
###########################

#first extract the barcodes
sel = d$seurat_clusters == 0 | d$seurat_clusters == 1 | d$seurat_clusters == 2 | d$seurat_clusters == 6 | d$seurat_clusters == 7
tcellBarcodes = colnames(gn10x)[sel]
sel = d$seurat_clusters == 3 | d$seurat_clusters == 5
monoBarcodes = colnames(gn10x)[sel]

#load the bug
tCounts = bug[(bug$barcode %in% tcellBarcodes) & (bug$gene == "ALDH2"),]$count
mCounts = bug[(bug$barcode %in% monoBarcodes) & (bug$gene == "ALDH2"),]$count
mCountsFilt = mCounts[mCounts < 30]
hend = max(c(mCountsFilt,tCounts))
numTTotMolecules = nrow(bug[bug$barcode %in% tcellBarcodes,])
numMTotMolecules = nrow(bug[bug$barcode %in% monoBarcodes,])
numTTotCounts = sum(bug[bug$barcode %in% tcellBarcodes,]$count)
numMTotCounts = sum(bug[bug$barcode %in% monoBarcodes,]$count)


ht = hist(tCounts, breaks=seq(0.5, hend+0.5, by=1), plot = FALSE)$counts
hm = hist(mCountsFilt, breaks=seq(0.5, hend+0.5, by=1), plot = FALSE)$counts

htScaled = ht/numTTotMolecules * 10^6
hmScaled = hm/numMTotMolecules * 10^6


dsPlot = data.frame(x = 1:hend, y = htScaled)
pC = ggplot(dsPlot,aes(x=x,y=y)) +
  geom_bar(stat="identity", fill = colors[1]) +
  labs(y="CPM", x="Counts per UMI", title="ALDH2, T cells") +
  ylim(0,50) +
  theme(panel.background = element_rect("white", "white", 0, 
                                        0, "white"),
        plot.title = element_text(face = "bold")
  )
pC

dsPlot = data.frame(x = 1:hend, y = hmScaled)
pD = ggplot(dsPlot,aes(x=x,y=y)) +
  geom_bar(stat="identity", fill = colors[2]) +
  labs(y="CPM", x="Counts per UMI", title="ALDH2, Monocytes") +
  ylim(0,50) +
  theme(panel.background = element_rect("white", "white", 0, 
                                        0, "white"),
        plot.title = element_text(face = "bold")
  )
pD


##############################
#Plot all variable genes
##############################

#first scale all genes to the same scale (sum of 1000)
scaled = TPMsForVariableGenes*1000/rowSums(TPMsForVariableGenes)
rowSums(scaled)
logScaled = log2(scaled+1)
cu = NULL
scLogExpr = NULL

for (i in 1:nrow(cuForVariableGenes)) {
  sel = !is.na(cuForVariableGenes[i,])
  cu = c(cu, as.numeric(cuForVariableGenes[i,sel]))  
  scLogExpr = as.numeric(c(scLogExpr, logScaled[i,sel]))
}

fit = lm(scLogExpr~cu)
summary(fit)# p-value: < 2.2e-16, F-Test

dsPlot = tibble(x=cu, y=scLogExpr)
dsText = tibble(x=1, y=10)


pG = ggplot2::ggplot(dsPlot,ggplot2::aes(x=x,y=y)) +
  ggplot2::geom_point(ggplot2::aes(x=x,y=y), color="black", shape=1, size=1) +
  geom_abline(slope = fit$coefficients[2], intercept = fit$coefficients[1], colour="#008800", size=1.3) +
  geom_text(data=dsText, label=paste0("R = ",format(cor(cu,scLogExpr), digits=2)), size=5, hjust = 0, parse=FALSE) +
  ggplot2::labs(y=expression(Log[2]*"(Norm. Expr + pc)"), x="Copies per UMI", title="Ampl. acr.clusters, uncorr.") +
  theme(panel.background = element_rect("white", "white", 0, 0, "white"),
        legend.position= "bottom", legend.direction = "horizontal",#, legend.title = element_blank())
        strip.text.x = element_text(size = 12, face = "bold"),
        #legend.position= "none",
        plot.title = element_text(face = "bold"),
        strip.background = element_blank())
pG


######################
# Now predict and see if the correlation goes down
######################

dsid = "EVALPBMC"
bug = getBug(dsid)

#first filter on the genes in the plot to make the bug smaller
genes = rownames(cuForVariableGenes)
bugFilt1 = bug[bug$gene %in% genes,]

predDS = UMIsForVariableGenes #just allocate the right size
predDS[,] = NA
predZTNB = UMIsForVariableGenes #just allocate the right size
predZTNB[,] = NA


for (i in 1:12) {
  print(i)
  sel = clust == i-1
  cellIds = names(clust)[sel]
  clustBug = bugFilt1[bugFilt1$barcode %in% cellIds,] #filter a bit more
  if (cellsPerCluster[i] > 300) {
    for (g in 1:length(genes)) {
      if (UMIsForVariableGenes[g,i] >= 20) {
        #Get the histogram
        counts = clustBug[clustBug$gene == genes[g],]$count
        h = hist(counts, breaks=seq(0.5, max(counts)+0.5, by=1), plot = FALSE)
        freq = h$mids
        counts = h$counts
        added = 0
        #preseq cannot handle if we have only ones, so modify the histogram slightly
        if ((length(freq)==1) & (freq[1] == 1)) {
          added = 2
          freq = c(1,2)
          counts = c(counts[1]+1,1)#room for improvement here
        }
        dd = as.matrix(data.frame(freq,counts));
        rSACZTNB = mod.ztnb.rSAC(dd, incTol = 1e-5, iterIncTol = 200);
        rSACDS = ds.rSAC(dd,mt=2)
        newCountsZTNB = rSACZTNB(10^20)
        newCountsDS = rSACDS(10^20)
        newCountsZTNB[newCountsZTNB < 0] = 0
        newCountsDS[newCountsDS < 0] = 0
        predZTNB[g,i] = newCountsZTNB - added;
        predDS[g,i] = newCountsDS - added;
      }
    }
  }
}


#remove all columns with just NA
predFiltZTNB = predZTNB[,colSums(!is.na(predZTNB)) != 0]
predLogZTNB = predFiltZTNB;
predFiltDS = predDS[,colSums(!is.na(predDS)) != 0]
predLogDS = predFiltDS;
#scale to the same number of average counts per cluster as before
colSumsData = TPMsForVariableGenes
colSumsData[is.na(predFiltZTNB)] = NA
scaleSizeZTNB = mean(colSums(colSumsData,na.rm=TRUE))
colSumsData = TPMsForVariableGenes
colSumsData[is.na(predFiltDS)] = NA
scaleSizeDS = mean(colSums(colSumsData,na.rm=TRUE))


for (i in 1:ncol(predFiltZTNB)) {
  TPMish = predFiltZTNB[,i]*scaleSizeZTNB/sum(predFiltZTNB, na.rm = TRUE)
  predLogZTNB[,i] = log2(TPMish + 1)
  TPMish = predFiltDS[,i]*scaleSizeDS/sum(predFiltDS, na.rm = TRUE)
  predLogDS[,i] = log2(TPMish + 1)
}

#use DS in main figure

#Get the average TPM per cluster before prediction
sumTPMBefore = sum(TPMsForVariableGenes[rownames(TPMsForVariableGenes) == "ALDH2",])

geneSel = rownames(predFiltDS) == "ALDH2" #ALDH2
cu = cuForVariableGenes[geneSel,]
pseudoTPM = predFiltDS[geneSel,] * sumTPMBefore/sum(predFiltDS[geneSel,])
logExpr = log2(pseudoTPM + 1)
fit = lm(logExpr~cu)
summary(fit) # p value: 0.5355, F-Test

dsPlot = tibble(x=cu, y=logExpr, col=c(colors[1],colors[1],colors[1],colors[2],"#000000",colors[2],colors[1],colors[1]))
dsText = tibble(x=cu + 0.2, y=logExpr)
#fix the texts a bit so they are visible and don't overlap
dsText$x[6] = dsText$x[6] - 1.5
dsText$y[1] = dsText$y[1] - 0.1

pF = ggplot2::ggplot(dsPlot,ggplot2::aes(x=x,y=y)) +
  ggplot2::geom_point(ggplot2::aes(x=x,y=y), color=dsPlot$col, shape=1, size=1) +
  geom_abline(slope = fit$coefficients[2], intercept = fit$coefficients[1], colour="#008800", size=1.3) +
  geom_text(data=dsText, label=texts, size=4, hjust = 0, color=dsPlot$col, parse=FALSE) +
  ggplot2::labs(y=expression(Log[2]*"(pseudoTPM + pc)"), x="Copies per UMI", title="ALDH2 acr. clusters, corr.") +
  theme(panel.background = element_rect("white", "white", 0, 0, "white"),
        legend.position= "bottom", legend.direction = "horizontal",#, legend.title = element_blank())
        strip.text.x = element_text(size = 12, face = "bold"),
        plot.title = element_text(face = "bold"),
        #legend.position= "none",
        strip.background = element_blank())
pF




#test for all variable genes

#first scale all genes to the same scale (sum of 1000)
scaled = predFiltDS*1000/rowSums(predFiltDS, na.rm=TRUE)
rowSums(scaled, na.rm=TRUE)
logScaled = log2(scaled+1)
cu = NULL
scLogExpr = NULL

for (i in 1:nrow(cuForVariableGenes)) {
  sel = !is.na(cuForVariableGenes[i,])
  cu = c(cu, as.numeric(cuForVariableGenes[i,sel]))  
  scLogExpr = as.numeric(c(scLogExpr, logScaled[i,sel]))
}

fit = lm(scLogExpr~cu)
summary(fit)#p = 5.111e-09, F-Test
dsPlot = tibble(x=cu, y=scLogExpr)
dsText = tibble(x=1, y=10)


pH = ggplot2::ggplot(dsPlot,ggplot2::aes(x=x,y=y)) +
  ggplot2::geom_point(ggplot2::aes(x=x,y=y), color="black", shape=1, size=1) +
  geom_abline(slope = fit$coefficients[2], intercept = fit$coefficients[1], colour="#008800", size=1.3) +
  geom_text(data=dsText, label=paste0("R = ",format(cor(cu,scLogExpr), digits=2)), size=5, hjust = 0, parse=FALSE) +
  ggplot2::labs(y=expression(Log[2]*"(Norm. Expr + pc)"), x="Copies per UMI", title="Ampl. acr. clusters, corr.") +
  theme(panel.background = element_rect("white", "white", 0, 0, "white"),
        legend.position= "bottom", legend.direction = "horizontal",#, legend.title = element_blank())
        strip.text.x = element_text(size = 12, face = "bold"),
        plot.title = element_text(face = "bold"),
        #legend.position= "none",
        strip.background = element_blank())
pH

fig5CH = ggarrange(pC, pD, pE, pF, pG, pH, nrow=2, ncol=3, labels=c("C","D","E","F","G","H"))

ggsave(
  paste0(figure_path, "Fig5C-H.png"),
  plot = fig5CH, device = "png",
  width = 9, height = 6, dpi = 300)

#Create supplementary figure with ZTNB
################

scaled = predFiltZTNB*1000/rowSums(predFiltZTNB, na.rm=TRUE)
rowSums(scaled, na.rm=TRUE)
logScaled = log2(scaled+1)
cu = NULL
scLogExpr = NULL

for (i in 1:nrow(cuForVariableGenes)) {
  sel = !is.na(cuForVariableGenes[i,])
  cu = c(cu, as.numeric(cuForVariableGenes[i,sel]))  
  scLogExpr = as.numeric(c(scLogExpr, logScaled[i,sel]))
}

fit = lm(scLogExpr~cu)
summary(fit)#p=0.4397, F-Test
dsPlot = tibble(x=cu, y=scLogExpr)
dsText = tibble(x=1, y=10)


pSup2 = ggplot2::ggplot(dsPlot,ggplot2::aes(x=x,y=y)) +
  ggplot2::geom_point(ggplot2::aes(x=x,y=y), color="black", shape=1, size=1) +
  geom_abline(slope = fit$coefficients[2], intercept = fit$coefficients[1], colour="#008800", size=1.3) +
  geom_text(data=dsText, label=paste0("R = ",format(cor(cu,scLogExpr), digits=2)), size=5, hjust = 0, parse=FALSE) +
  ggplot2::labs(y=expression(Log[2]*"(Norm. Expr + pc)"), x="Copies per UMI") +
  theme(panel.background = element_rect("white", "white", 0, 0, "white"),
        legend.position= "bottom", legend.direction = "horizontal",#, legend.title = element_blank())
        strip.text.x = element_text(size = 12, face = "bold"),
        #legend.position= "none",
        strip.background = element_blank())
pSup2

ggsave(
  paste0(figure_path, "FigS_AcrossClustersZTNB_5C-H.png"),
  plot = pSup2, device = "png",
  width = 3, height = 3, dpi = 300)

