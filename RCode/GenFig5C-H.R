#before anything else, you need to setup paths:
#for example
#source("C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/paths.R")
#or, in colab:
#source("GRNP_2020/RCode/pathsGoogleColab.R")

source(paste0(sourcePath,"ButterflyHelpers.R"))
source(paste0(sourcePath,"CCCHelpers.R"))
source(paste0(sourcePath,"preseqHelpers.R"))
#source(paste0(sourcePath,"Sm2Helpers.R"))
sourcePath2 = "C:/Work/R/ButterflyPublicationRepo/GRNP_2020/RCode/"
source(paste0(sourcePath2,"Sm2Helpers.R"))
source(paste0(sourcePath,"ggplotHelpers.R"))

EVALPBMC_SM2_path = "D:/EVALPBMC-SM2/output/"

library(ggplot2)
library("ggpubr")
library(tidyverse)
library(BUSpaRse)
library("matrixStats")

colors = gg_color_hue(2)


########################
# Load data
########################

#Get tx data from EVALPBMC 10X
tx10x = read_count_output(dir="E:/Butterfly/EVALPBMC/cu_processing/count_output_tx_m", "output", tcc=TRUE)

dim(tx10x)
countsPerCell = sparse_Sums(tx10x)
tx10x = tx10x[,countsPerCell >= 200]
tx10x = tx10x[sparse_Sums(tx10x, rowSums = TRUE) != 0,]
dim(tx10x)#127593   5176   _m:2590285    5319


#Seurat cannot handle the tx data
#So, cluster on gene data instead, it doesn't really matter
gn10x = read_count_output(dir="E:/Butterfly/EVALPBMC/cu_processing/count_output", "output", tcc=FALSE)
countsPerCell = sparse_Sums(gn10x)
gn10x = gn10x[,countsPerCell >= 200]
gn10x = gn10x[sparse_Sums(gn10x, rowSums = TRUE) != 0,]
gn10xRawCounts = read_count_output(dir="E:/Butterfly/EVALPBMC/cu_processing/raw_count_output/", "output", tcc=FALSE)
gn10xRawCounts = gn10xRawCounts[,countsPerCell >= 200]
gn10xRawCounts = gn10xRawCounts[sparse_Sums(gn10xRawCounts, rowSums = TRUE) != 0,]


set.seed(1)

library(Seurat)

#run through Seurat
d = CreateSeuratObject(counts = gn10x, project = "10x", min.cells = 0, min.features = 0)

#create variable for mitochondrial count fraction
d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "^MT-")

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

d = RenameIdents(object = d, `0` = "T1", `1` = "T2", `2` = "T3", `3` = "M1", `4` = "B", `5` = "T4", `6` = "M2", `7` = "U1", `8` = "U2", `9` = "U3", `10` = "U4",  `11` = "U5")

#show clustering
pSupA = DimPlot(d, reduction = "umap")


#FeaturePlot(d, features = c("ENST00000416293.7","ENST00000261733.6","ENST00000548536.1","ENST00000549106.1","ENST00000551450.1"))

ggsave(
  paste0(figure_path, "FigS_Seurat.png"),
  plot = pSupA, device = "png",
  width = 4, height = 4, dpi = 300)




#library(DescTools)
library(textTinyR)

#calculate counts per UMI
cuMatrix = gn10x

for (i in 1:ncol(cuMatrix)) {
  if (i %% 100 == 0) print(i)
  sel = cuMatrix[,i] != 0
  cuMatrix[sel,i] = gn10xRawCounts[sel,i]/gn10x[sel,i]
}


######################
# Compare CU across clusters
######################

cuPerCluster = matrix(, nrow = nrow(cuMatrix), ncol = 12) #number of clusters is 12
rownames(cuPerCluster) = rownames(gn10x)
UMIsPerCluster = matrix(, nrow = nrow(cuMatrix), ncol = 12) #number of clusters is 12
rownames(UMIsPerCluster) = rownames(gn10x)
cellsPerCluster = rep(NA,12)
clusterNames = c("T1","T2","T3","M1","B","T4","M2","7","8","9","10","11") #cluster 4 is b cells, 3,6 are monocytes, 0,1,2,5 are T cells/NK cells
colnames(cuPerCluster) = clusterNames
colnames(UMIsPerCluster) = clusterNames

for (i in 1:12) {
  sel = d$seurat_clusters == i-1
  cellsPerCluster[i] = sum(sel)
  cuPerCluster[,i] = sparse_Sums(gn10xRawCounts[,sel], rowSums = TRUE)/sparse_Sums(gn10x[,sel], rowSums = TRUE)
  UMIsPerCluster[,i] = sparse_Sums(gn10x[,sel], rowSums = TRUE)
}
cellsPerCluster
clusterFilt = 0:11
clusterFilt = clusterFilt[cellsPerCluster > 100]

cuPerClusterFilt = cuPerCluster[,cellsPerCluster > 300]
UMIsPerClusterFilt = UMIsPerCluster[,cellsPerCluster > 300]
cuPerClusterFilt[UMIsPerClusterFilt < 20] = NA

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

geneSel = rownames(cuForVariableGenes) == "ENSG00000111275.12" #ALDH2
cu = cuForVariableGenes[geneSel,]
logExpr = logtrans[geneSel,]
texts = names(cu)

dsPlot = tibble(x=cu, y=logExpr, col=c(colors[1],colors[1],colors[1],colors[2],"#000000",colors[1],colors[2]))
dsText = tibble(x=cu + 0.2, y=logExpr)
fit = lm(logExpr~cu)
#fix the texts a bit so they are visible and don't overlap
dsText$x[7] = dsText$x[7] - 1.5
dsText$y[2] = dsText$y[2] - 0.05
dsText$y[6] = dsText$y[6] + 0.05

pA = ggplot2::ggplot(dsPlot,ggplot2::aes(x=x,y=y)) +
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
pA


###########################
#look at the histograms for T cells and Monocytes
###########################

#first extract the barcodes
sel = d$seurat_clusters == 0 | d$seurat_clusters == 1 | d$seurat_clusters == 2 | d$seurat_clusters == 5
tcellBarcodes = colnames(gn10x)[sel]
sel = d$seurat_clusters == 3 | d$seurat_clusters == 6
monoBarcodes = colnames(gn10x)[sel]

#load the bug
loadBug("EVALPBMC")
bug = getBug("EVALPBMC")
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
pB = ggplot(dsPlot,aes(x=x,y=y)) +
  geom_bar(stat="identity", fill = colors[1]) +
  labs(y="CPM", x="Counts per UMI", title="ALDH2, T cells") +
  ylim(0,50) +
  theme(panel.background = element_rect("white", "white", 0, 
                                        0, "white"),
        plot.title = element_text(face = "bold")
  )
pB

dsPlot = data.frame(x = 1:hend, y = hmScaled)
pC = ggplot(dsPlot,aes(x=x,y=y)) +
  geom_bar(stat="identity", fill = colors[2]) +
  labs(y="CPM", x="Counts per UMI", title="ALDH2, Monocytes") +
  ylim(0,50) +
  theme(panel.background = element_rect("white", "white", 0, 
                                        0, "white"),
        plot.title = element_text(face = "bold")
  )
pC

#calculate CPM of raw counts
tRawCountsTPM = sum(tCounts) * 10^6/numTTotCounts #7.13283
mRawCountsTPM = sum(mCounts) * 10^6/numMTotCounts #196.5003

############################
## Now investigate ALDH2 in Smart-Seq2 data
############################
scEvSm2List = ReadKallistoSm2TxMatrix(paste0(dataPath, "EVALPBMC/bus_output/transcripts_to_genes.txt"),
                                      "E:/Butterfly/EVALPBMC/bulk/PBMC1/abundance.tsv", EVALPBMC_SM2_path)
countsMatrix = scEvSm2List[[1]]
tpmMatrix = scEvSm2List[[2]]

#filter out zero genes
countsMatrix = countsMatrix[rowSums(countsMatrix) != 0,]
tpmMatrix = tpmMatrix[rowSums(tpmMatrix) != 0,]
dim(countsMatrix)

#run through Seurat
library(Seurat)
set.seed(1)

d = CreateSeuratObject(counts = tpmMatrix, project = "sm2", min.cells = 0, min.features = 0)

#create variable for mitochondrial count fraction
#d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "^MT-")

#library size normalization + change to log scale
d <- NormalizeData(d, normalization.method = "LogNormalize", scale.factor = 10000)

#Finds the most variable genes (takes a few minutes to run)
d <- FindVariableFeatures(d, selection.method = "vst", nfeatures = 2000)


#Make mean and variance the same for all genes:
d <- ScaleData(d, features = rownames(d))

#Principal component analysis
d <- RunPCA(d, features = VariableFeatures(object = d))

#do clustering
d <- FindNeighbors(d, dims = 1:10)
d <- FindClusters(d, resolution = 0.5)

#Generate UMAP map (similar to t-SNE, but usually better)
d <- RunUMAP(d, dims = 1:10)

#show clustering
DimPlot(d, reduction = "umap")

sm2Clust = d$seurat_clusters

#figure out what the clusters are
#b cells
#variants of CD19:
#ENST00000324662.7	ENSG00000177455.12	CD19
#ENST00000565089.5	ENSG00000177455.12	CD19
#ENST00000567541.5	ENSG00000177455.12	CD19
#ENST00000566890.1	ENSG00000177455.12	CD19
#ENST00000567368.1	ENSG00000177455.12	CD19
#ENST00000538922.5	ENSG00000177455.12	CD19
FeaturePlot(d, features = c("ENST00000324662.7", "ENST00000565089.5", "ENST00000567541.5", "ENST00000566890.1", "ENST00000567368.1", "ENST00000538922.5")) #_, CD19,
#so, cluster 3 is B cells

#monocytes
#ENST00000261267.6	ENSG00000090382.6	LYZ
#ENST00000549690.1	ENSG00000090382.6	LYZ
#ENST00000548839.1	ENSG00000090382.6	LYZ

FeaturePlot(d, features = c("ENST00000261267.6", "ENST00000549690.1", "ENST00000548839.1"))
#cluster 1 is monocytes


#T cells/NK cells
#ENST00000300692.8	ENSG00000167286.9	CD3D
#ENST00000534687.5	ENSG00000167286.9	CD3D
#ENST00000529594.5	ENSG00000167286.9	CD3D
#ENST00000392884.2	ENSG00000167286.9	CD3D
#ENST00000526561.1	ENSG00000167286.9	CD3D

FeaturePlot(d, features = c("ENST00000300692.8", "ENST00000534687.5", "ENST00000529594.5", "ENST00000392884.2", "ENST00000526561.1"))
#cluster 0 is t cells


#ALDH2
#ENST00000416293.7	ENSG00000111275.12	ALDH2
#ENST00000261733.6	ENSG00000111275.12	ALDH2
#ENST00000548536.1	ENSG00000111275.12	ALDH2
#ENST00000549106.1	ENSG00000111275.12	ALDH2
#ENST00000551450.1	ENSG00000111275.12	ALDH2
FeaturePlot(d, features = c("ENST00000416293.7","ENST00000261733.6","ENST00000548536.1","ENST00000549106.1","ENST00000551450.1"))

#ok, so we are interested in comparing the T cells and monocytes
monoPooled = rowMeans(tpmMatrix[,sm2Clust == 1])
tPooled = rowMeans(tpmMatrix[,sm2Clust == 0])

#TPM them again just to be safe
monoPooled = monoPooled*10^6/sum(monoPooled)
tPooled = tPooled*10^6/sum(tPooled)

ALDH2Transcr = c("ENST00000416293.7","ENST00000261733.6","ENST00000548536.1","ENST00000549106.1","ENST00000551450.1")

monoPooledFilt = monoPooled[names(monoPooled) %in% ALDH2Transcr]
tPooledFilt = tPooled[names(tPooled) %in% ALDH2Transcr]
#The gene expression matches reasonably between smart-seq2 and raw counts
sum(monoPooledFilt) #298.4646
mRawCountsTPM #196.5003
sum(tPooledFilt) #5.049317
tRawCountsTPM #7.13283

#ENST00000551450.1 seems to be a good candidate for the lowly copied transcript. Can we estimate this using the 10X bug?



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

dsPlot = tibble(x=cu, y=scLogExpr)
dsText = tibble(x=1, y=10)


pE = ggplot2::ggplot(dsPlot,ggplot2::aes(x=x,y=y)) +
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
pE


######################
# Now predict and see if the correlation goes down
######################

dsid = "EVALPBMC"
loadBug(dsid, "ENS")
bugEns = getBug(dsid,"ENS")

#first filter on the genes in the plot to make the bug smaller
genes = rownames(cuForVariableGenes)
bugFilt1 = bugEns[bugEns$gene %in% genes,]

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
colSumsData[is.na(predFilt)] = NA
scaleSize = mean(colSums(colSumsData,na.rm=TRUE))


for (i in 1:ncol(predFiltZTNB)) {
  TPMish = predFiltZTNB[,i]*scaleSize/sum(predFiltZTNB, na.rm = TRUE)
  predLogZTNB[,i] = log2(TPMish + 1)
  TPMish = predFiltDS[,i]*scaleSize/sum(predFiltDS, na.rm = TRUE)
  predLogDS[,i] = log2(TPMish + 1)
}

#use DS in main figure

#Get the average TPM per cluster before prediction
sumTPMBefore = sum(TPMsForVariableGenes[rownames(TPMsForVariableGenes) == "ENSG00000111275.12",])

geneSel = rownames(predFiltDS) == "ENSG00000111275.12" #ALDH2
cu = cuForVariableGenes[geneSel,]
pseudoTPM = predFiltDS[geneSel,] * sumTPMBefore/sum(predFiltDS[geneSel,])
logExpr = log2(pseudoTPM + 1)
fit = lm(logExpr~cu)

dsPlot = tibble(x=cu, y=logExpr, col=c(colors[1],colors[1],colors[1],colors[2],"#000000",colors[1],colors[2]))
dsText = tibble(x=cu + 0.2, y=logExpr)
#fix the texts a bit so they are visible and don't overlap
dsText$x[7] = dsText$x[7] - 1.5
dsText$y[2] = dsText$y[2] - 0.05
dsText$y[3] = dsText$y[3] + 0.05

pD = ggplot2::ggplot(dsPlot,ggplot2::aes(x=x,y=y)) +
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
pD




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
dsPlot = tibble(x=cu, y=scLogExpr)
dsText = tibble(x=1, y=10)


pF = ggplot2::ggplot(dsPlot,ggplot2::aes(x=x,y=y)) +
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
pF

fig5CH = ggarrange(pB, pC, pA, pD, pE, pF, nrow=2, ncol=3, labels=c("C","D","E","F","G","H"))

ggsave(
  paste0(figure_path, "Fig5C-H.png"),
  plot = fig5CH, device = "png",
  width = 9, height = 6, dpi = 300)

#Create potential supplementary figure with ZTNB
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
  paste0(figure_path, "FigS_AcrossClustersZTNB.png"),
  plot = pSup2, device = "png",
  width = 3, height = 3, dpi = 300)




########################################
# Investigate splice variants
########################################

###################
# 10x - see if we can find any differences in splice variants across clusters
###################

#load the bus file (Not bug)
bus = read.table(paste0(dataPath,"EVALPBMC/cu_processing/bus.txt"), stringsAsFactors = F)

#load transcripts file
trs = read.table(paste0(dataPath,"EVALPBMC/cu_processing/transcripts.txt"), stringsAsFactors = F)
#get transcript ids (0-based indices in the transcripts file)
ALDH2Transcr = c("ENST00000416293.7","ENST00000261733.6","ENST00000548536.1","ENST00000549106.1","ENST00000551450.1")
trIds = which(trs[,1] %in% ALDH2Transcr) -1 #10904 10905 10906 10907 10908

ecMatrix = read.table(paste0(dataPath,"EVALPBMC/cu_processing/matrix.ec"), stringsAsFactors = F)
#get all ec:s connected to any of these
#first filter all that doesn't contain any of those strings, to reduce the workload
sel = rep(FALSE, nrow(ecMatrix))
for (id in trIds) {
  sel = sel | grepl(as.character(id), ecMatrix$V2, fixed = TRUE)  
}
sum(sel)
ecFilt1 = ecMatrix[sel,]
dim(ecFilt1) #2773    2 - some of these are not real matches, but could be 107 in 1107 and similar
#create a dataset with ec:s and transcripts, where all other transcripts than the ALDH2 have been removed
expandedEC = separate_rows(ecFilt1, V2, convert = T)
dim(expandedEC)#305626      2
#Now remove all rows that doesn't contain ALDH2 transcripts
expandedECFilt = expandedEC[expandedEC$V2 %in% trIds,]
dim(expandedECFilt)#2350    2

collapsedEC = expandedECFilt %>% group_by(V1) %>% do(transcriptsList=c(.$V2))
dim(collapsedEC)
colnames(collapsedEC) = c("x","y")

#Ok, the trick is now to only keep the BUS records that have the same UMI and barcode as BUG records that have ALDH2 as gene.
#we know that these are uniquely mapped

aldh2Bug = bug[bug$gene == "ALDH2",]
dim(aldh2Bug)#851   4
aldhKeys = paste0(aldh2Bug$barcode, aldh2Bug$UMI)

#first filter to get a smaller dataset to work with
busFilt1 = bus[bus$V2 %in% aldh2Bug$UMI,]

busKeys = paste0(busFilt1$V1, busFilt1$V2)

busFilt2 = busFilt1[busKeys %in% aldhKeys,]
dim(busFilt2) #1186    4

keyedBus = tibble(key=paste(busFilt2$V1,busFilt2$V2), transcripts = lookup(busFilt2$V3, collapsedEC), count = busFilt2$V4 )

transcriptIntercept = function(transcripts) {
  intersection = transcripts[1][[1]]
  if (length(transcripts) > 1) {
    for (i in 2:length(transcripts)) {
      intersection = intersect(intersection, transcripts[i][[1]])
    }
  }
  return(list(intersection))
}

transcriptIntercept(list(c("A", "B"), c("B")))
transcriptIntercept(list(c("A", "B")))

intersect(c("A", "B"), c("B"))

xx  = list(c("A", "B"))

xx[1][[1]]

transcriptIntercept(list(c("A", "B")))

testBus = keyedBus[1:10,]
pooledBus = testBus %>% group_by(key) %>% summarize(transcripts2 = transcriptIntercept(transcripts), count2 = sum(count))


pooledBus = keyedBus %>% group_by(key) %>% summarize(transcripts2 = transcriptIntercept(transcripts), count2 = sum(count))
dim(pooledBus)#851   3, good, same as from bug file
numTr = sapply(pooledBus$transcripts2, FUN=length) #almost only unique reads - filter the others
pooledBusFilt = pooledBus[numTr == 1,]
dim(pooledBusFilt)#840   3
#now unlist the transcript
unlist(pooledBusFilt$transcripts2)
sum(pooledBusFilt$count2 == 1)/length(pooledBusFilt$count2)










#hmm, not sure how to do this... we have ec:s, not transcripts....
geneSel = rownames(cuForVariableGenes) == "ENSG00000111275.12" #ALDH2
cuForVariableGenes[geneSel,] #cluster 3 and 6 (monocytes/DC) have high CU, T cells low, B cells not quite as low

#transcript bustools id
#ENST00000416293.7 - 10904
#ENST00000261733.6 - 10905
#ENST00000548536.1 - 10906
#ENST00000549106.1 - 10907
#ENST00000551450.1 - 10908

#ECs referring to those transcripts
#EC     Transcript
#10904	10904
#10905	10905
#10906	10906
#10907	10907
#10908	10908
#192473	10904,10905
#210811	10904,10905,10906
#227292	10904,10905,10906,10909
#548432	10904,10905,10909
#677573	10904,10905,10906,10907
#699527	10904,10905,10907
#191329	1232,1233,1234,1235,1753,2126,2387,2607,2608,2697,4513,4515,4669,4906,6910,6919,7065,7066,7080,7085,7203,7557,7714,7773,7774,8197,8202,8638,8639,8709,9291,9331,9362,9436,9596,9800,9987,10577,10905,11049,11308,11325,11335,11764,13246,13378,13440,14547,14548,14554,15111,15112,15182,15183,15184,15491,15709,15721,15769,16266,16371,19330,19549,20203,20207,20496,21046,21047,21067,21916,22233,22400,22545,22546,23797,23803,23970,24316,25410,25550,25744,25908,26084,26788,26933,26934,26935,27040,27464,28338,28980,29357,29358,29580,30853,31271,31319,31453,31458,31694,31716,32066,32319,32378,32967,32968,33309,33801,34239,34886,34947,35602,35775,36417,37436,38065,38648,38654,38656,39230,39418,39419,40319,40357,40850,40851,41405,41731,41827,41828,42756,43181,43730,45613,46986,49102,50352,50792,51496,51555,52699,53648,53943,55194,55907,55915,56131,56132,56600,57474,57476,57792,57823,59048,59161,59294,59781,59867,60716,60718,61066,62187,62231,62351,62892,63432,63450,63457,63765,63983,63985,64113,64751,65134,65385,65387,66393,67130,67337,67708,67810,67881,68718,69243,69745,69790,69899,69903,69904,69912,69914,70026,70298,70454,70593,70679,70680,70873,70874,70923,71204,71290,71414,71905,72047,73284,73430,73612,73924,75203,75483,76084,77202,77220,77256,77412,77438,77550,78086,78680,78681,79267,80210,80317,80848,81346,81443,81455,81530,81736,82614,83166,83193,83865,83935,84580,84581,84582,84583,85317,85420,85421,85637,85817,85818,85998,86173,86667,86669,86670,88512,88874,88995,89489,89676,89716,89717,90214,90220,91068,91492,91498,92049,92243,92274,92668,92996,93113,93114,93199,93200,93727,94136,94579,95636,95924,96115,96510,96929,96930,97179,97180,97313,97401,97808,97811,97830,97852,97920,97931,98027,98758,100403,100456,101096,101311,101367,101398,101400,101406,102689,104052,106217,106222,109895,109907,109908,110106,110230,110587,110794,111434,111580,111581,111682,112619,113559,114113,114977,114978,115145,115931,115935,116044,116317,116318,116847,117416,117935,118050,118427,119674,119678,119876,119880,122386,122704,122907,123284,123288,123708,123754,123756,123759,123892,123893,124265,124266,124267,124434,125131,125169,125619,125635,126159,126719,127670,128501,128687,128762,128843,129562,129842,130003,130004,130145,131183,131188,131290,131999,132000,132001,132007,132008,132009,132655,133345,133638,133644,133720,134564,136878,137806,138160,138251,138791,138796,138797,138798,138799,138800,138801,138802,138803,138904,138905,139056,139057,139283,139892,140299,140872,140873,141193,141207,141217,141218,141241,141403,141404,141505,141508,142345,142350,144113,146226,146228,146512,146571,147268,147420,147457,147564,147565,147566,148099,148513,148591,149048,149410,149426,149848,149949,150390,150493,151190,151210,151447,151776,151985,152059,152122,152337,152596,152895,153194,153312,153353,153555,154027,155157,155714,155756,156053,156119,156397,156399,156609,156610,157571,158139,159323,159330,159718,159727,160050,160051,160208,160814,161389,161642,161906,161907,161917,162562,162591,163377,164135,164189,164190,164316,165275,165276,165510,165511,166080,166081,166201,166471,166473,167205,168185,169435,169438,170121,170175,170199,170355,171213,171486,172584,172681,176215,179237,179965,181223,183242,183368,186147,187145
#221334	10905,13733,87956,111730,113166,128816,130723,157999
#230158	10905,12109
#243407	10905,10906,10909
#253658	10905,126899
#260074	1170,1174,10905,26558,26664,26667,52353,55179,55180,63670,72860,78215,81261,81262,81264,91774,103384,103387,126967,133016
#265257	10905,25591,25594
#270467	10905,40964


#c("ENST00000416293.7","ENST00000261733.6","ENST00000548536.1","ENST00000549106.1","ENST00000551450.1")
aldh2_10X_tx = tx10x[rownames(tx10x) %in% c(10904,10905,10906,10907,10908),]
sparse_Sums(aldh2_10X_tx, rowSums=TRUE)


##################
#Smart-Seq2
##################

scEvSm2List = ReadKallistoSm2TxMatrix(paste0(dataPath, "EVALPBMC/bus_output/transcripts_to_genes.txt"),
                                      "E:/Butterfly/EVALPBMC/bulk/PBMC1/abundance.tsv", EVALPBMC_SM2_path)
dim(scEvSm2List)

countsMatrix = scEvSm2List[[1]]
tpmMatrix = scEvSm2List[[2]]

#filter out zero genes
countsMatrix = countsMatrix[rowSums(countsMatrix) != 0,]
tpmMatrix = tpmMatrix[rowSums(tpmMatrix) != 0,]
dim(countsMatrix)

#run through Seurat
library(Seurat)
set.seed(1)

d = CreateSeuratObject(counts = tpmMatrix, project = "sm2", min.cells = 0, min.features = 0)

#create variable for mitochondrial count fraction
d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "^MT-")

#library size normalization + change to log scale
d <- NormalizeData(d, normalization.method = "LogNormalize", scale.factor = 10000)

#Finds the most variable genes (takes a few minutes to run)
d <- FindVariableFeatures(d, selection.method = "vst", nfeatures = 2000)


#Make mean and variance the same for all genes:
d <- ScaleData(d, features = rownames(d))

#Principal component analysis
d <- RunPCA(d, features = VariableFeatures(object = d))

#do clustering
d <- FindNeighbors(d, dims = 1:10)
d <- FindClusters(d, resolution = 0.5)

#Generate UMAP map (similar to t-SNE, but usually better)
d <- RunUMAP(d, dims = 1:10)

#show clustering
DimPlot(d, reduction = "umap")

#figure out what the clusters are
#b cells
#variants of CD19:
#ENST00000324662.7	ENSG00000177455.12	CD19
#ENST00000565089.5	ENSG00000177455.12	CD19
#ENST00000567541.5	ENSG00000177455.12	CD19
#ENST00000566890.1	ENSG00000177455.12	CD19
#ENST00000567368.1	ENSG00000177455.12	CD19
#ENST00000538922.5	ENSG00000177455.12	CD19
FeaturePlot(d, features = c("ENST00000324662.7", "ENST00000565089.5", "ENST00000567541.5", "ENST00000566890.1", "ENST00000567368.1", "ENST00000538922.5")) #_, CD19,

#monocytes
#ENST00000261267.6	ENSG00000090382.6	LYZ
#ENST00000549690.1	ENSG00000090382.6	LYZ
#ENST00000548839.1	ENSG00000090382.6	LYZ

FeaturePlot(d, features = c("ENST00000261267.6", "ENST00000549690.1", "ENST00000548839.1"))

#T cells/NK cells
#ENST00000352580.7	ENSG00000153563.15	CD8A
#ENST00000283635.7	ENSG00000153563.15	CD8A
#ENST00000409511.6	ENSG00000153563.15	CD8A
#ENST00000409781.1	ENSG00000153563.15	CD8A
FeaturePlot(d, features = c("ENST00000352580.7", "ENST00000283635.7", "ENST00000409511.6", "ENST00000409781.1"))



#ALDH2
#ENST00000416293.7	ENSG00000111275.12	ALDH2
#ENST00000261733.6	ENSG00000111275.12	ALDH2
#ENST00000548536.1	ENSG00000111275.12	ALDH2
#ENST00000549106.1	ENSG00000111275.12	ALDH2
#ENST00000551450.1	ENSG00000111275.12	ALDH2
FeaturePlot(d, features = c("ENST00000416293.7","ENST00000261733.6","ENST00000548536.1","ENST00000549106.1","ENST00000551450.1"))

cuForVariableGenes[genes == "ENSG00000167604.14",]

#ENST00000590828.5	ENSG00000167604.14	NFKBID
#ENST00000641389.1	ENSG00000167604.14	NFKBID
#ENST00000606253.5	ENSG00000167604.14	NFKBID
#ENST00000586361.5	ENSG00000167604.14	NFKBID
#ENST00000588497.2	ENSG00000167604.14	NFKBID
#ENST00000590094.6	ENSG00000167604.14	NFKBID
#ENST00000591730.5	ENSG00000167604.14	NFKBID
#ENST00000585925.6	ENSG00000167604.14	NFKBID
#ENST00000588039.5	ENSG00000167604.14	NFKBID
#ENST00000585544.1	ENSG00000167604.14	NFKBID
#ENST00000396901.5	ENSG00000167604.14	NFKBID
FeaturePlot(d, features = c("ENST00000590828.5","ENST00000641389.1","ENST00000606253.5","ENST00000586361.5","ENST00000588497.2","ENST00000590094.6",
                            "ENST00000591730.5","ENST00000585925.6","ENST00000588039.5","ENST00000585544.1","ENST00000396901.5"))

#ENST00000368732.5	ENSG00000143546.9	S100A8
#ENST00000368733.3	ENSG00000143546.9	S100A8
#ENST00000477801.1	ENSG00000143546.9	S100A8
cuForVariableGenes[genes == "ENSG00000143546.9",]
FeaturePlot(d, features = c("ENST00000368732.5","ENST00000368733.3","ENST00000477801.1"))


cuForVariableGenes[genes == "ENSG00000168685.14",]
cuForVariableGenes[genes == "ENSG00000099337.4",]
#ENST00000263372.4	ENSG00000099337.4	KCNK6
#ENST00000588137.1	ENSG00000099337.4	KCNK6
FeaturePlot(d, features = c("ENST00000263372.4","ENST00000588137.1"))
cuForVariableGenes[genes == "ENSG00000109321.10",]
cuForVariableGenes[genes == "ENSG00000112096.17",]

#ENST00000538183.6	ENSG00000112096.17	SOD2
#ENST00000367055.8	ENSG00000112096.17	SOD2
#ENST00000546260.5	ENSG00000112096.17	SOD2
#ENST00000367054.6	ENSG00000112096.17	SOD2
#ENST00000444946.6	ENSG00000112096.17	SOD2
#ENST00000337404.8	ENSG00000112096.17	SOD2

#ENST00000535459.5	ENSG00000112096.17	SOD2
#ENST00000541573.5	ENSG00000112096.17	SOD2
#ENST00000540491.1	ENSG00000112096.17	SOD2
#ENST00000545162.5	ENSG00000112096.17	SOD2
#ENST00000535561.5	ENSG00000112096.17	SOD2
#ENST00000537657.5	ENSG00000112096.17	SOD2
#ENST00000401980.3	ENSG00000112096.17	SOD2
#ENST00000452684.2	ENSG00000112096.17	SOD2

FeaturePlot(d, features = c("ENST00000538183.6","ENST00000367055.8","ENST00000546260.5","ENST00000367054.6","ENST00000444946.6","ENST00000337404.8",
                            "ENST00000535459.5","ENST00000541573.5","ENST00000540491.1","ENST00000545162.5","ENST00000535561.5","ENST00000537657.5","ENST00000401980.3","ENST00000452684.2"))


