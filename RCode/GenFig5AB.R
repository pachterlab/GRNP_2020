#before anything else, you need to setup paths:
#for example
#source("C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/paths.R")
#or, in colab:
#source("GRNP_2020/RCode/pathsGoogleColab.R")

source(paste0(sourcePath,"ButterflyHelpers.R"))

library(ggplot2)
library("ggpubr")

source(paste0(sourcePath,"BinomialDownSampling.R"))
source(paste0(sourcePath,"CCCHelpers.R"))


#remember that the genes file is missing in the predquant output folder. Just copy the genes file output by collapse into the folder and rename it.


library(BUSpaRse)
library(textTinyR)
library(Matrix.utils)
library(Seurat)

################################################################
## Use the PBMC_NG + and PBMC_NG_2 datasets and batch correct
################################################################

#remember that the genes file is missing in the predquant output folder. Just copy the genes file output by collapse into the folder and rename it.

#read NG data
#############
NGMatrix = read_count_output(paste0(dataPath, "pbmc5kNextGem/bus_output/genecounts"), "output", tcc=F)
#filter cells with less than 200 reads
cs = sparse_Sums(NGMatrix, rowSums = F)

NGMatrix = NGMatrix[,cs > 200]

#convert the genes to gene symbols
tr2g = read.table(paste0(dataPath, "pbmc5kNextGem/bus_output/transcripts_to_genes.txt"), stringsAsFactors = F)
lookupTable = tr2g[,2:3]
lookupTable= unique(lookupTable)

outGenes = lookup(row.names(NGMatrix), lookupTable)

NGMat = NGMatrix
row.names(NGMat) = outGenes
NGMat = aggregate.Matrix(NGMat, outGenes, fun='sum')

#Read NG_2 data
###############
NG_2Matrix = read_count_output(paste0(dataPath, "pbmc_NG2/bus_output/genecounts"), "output", tcc=F)
#filter cells with less than 200 reads
cs = sparse_Sums(NG_2Matrix, rowSums = F)

NG_2Matrix = NG_2Matrix[,cs > 200]

#convert the genes to gene symbols
tr2g = read.table(paste0(dataPath, "pbmc_NG2/bus_output/transcripts_to_genes.txt"), stringsAsFactors = F)
lookupTable = tr2g[,2:3]
lookupTable= unique(lookupTable)

outGenes = lookup(row.names(NG_2Matrix), lookupTable)

NG_2Mat = NG_2Matrix
row.names(NG_2Mat) = outGenes
NG_2Mat = aggregate.Matrix(NG_2Mat, outGenes, fun='sum')

#Join the two datasets and process in Seurat
############################################
#now, test to merge (inner_join) the two sparse matrices
mergedmat = merge(NGMat, NG_2Mat, row.names(NGMat), row.names(NG_2Mat), all.x = F, all.y = F)

#verify that the merge seemed to work
sum(mergedmat[row.names(mergedmat) == "A2M",] == c(NGMat[row.names(NGMat) == "A2M",],NG_2Mat[row.names(NG_2Mat) == "A2M",]))
dim(mergedmat)#seems to have worked for that gene at least

#now, investigate the matrix in Seurat
#########################################

set.seed(1)

seurObj = CreateSeuratObject(counts = mergedmat, project = "merged", min.cells = 0, min.features = 0)
seurObj[["ds_source"]] = factor(as.numeric(c(rep(0, dim(NGMat)[2]), rep(1, dim(NG_2Mat)[2]))), c(0,1), c("PBMC_NG","PBMC_NG_2"))
seurObj[["percent.mt"]] <- PercentageFeatureSet(seurObj, pattern = "^MT-")
VlnPlot(seurObj, features = c("percent.mt"))

cellFilter = seurObj$percent.mt < 15
sourceFilt1 = seurObj$ds_source[cellFilter]
mergedmatFilt1 = mergedmat[,cellFilter] #cells already filtered at counts >= 200

#easiest to recreate the Seurat object:
seurObj = CreateSeuratObject(counts = mergedmatFilt1, project = "merged", min.cells = 0, min.features = 0)
seurObj[["ds_source"]] = sourceFilt1
seurObj <- NormalizeData(seurObj, normalization.method = "LogNormalize", scale.factor = 10000)
seurObj <- FindVariableFeatures(seurObj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurObj)
seurObj <- ScaleData(seurObj, features = all.genes)
seurObj <- RunPCA(seurObj, features = VariableFeatures(object = seurObj))
seurObj <- FindNeighbors(seurObj, dims = 1:10)
seurObj <- FindClusters(seurObj, resolution = 0.5)
seurObj <- RunUMAP(seurObj, dims = 1:10)
#DimPlot(seurObj, reduction = "umap", group.by = "ds_source", cols=c("blue", "red"))

FeaturePlot(seurObj, features = c("CD19", "CD8A", "CD3D")) # The clusters are T cells

supFigA = DimPlot(seurObj, reduction = "umap")
supFigA

ggsave(
  paste0(figure_path, "FigS_batch.png"),
  plot = supFigA, device = "png",
  width = 8, height = 5, dpi = 300)

DimPlot(seurObj, reduction = "umap", group.by = "ds_source", cols=c("blue", "red"))

#seurObj2 = subset(seurObj, subset = seurat_clusters == 0 | seurat_clusters == 3 | seurat_clusters == 5 | seurat_clusters == 6) 
#seurObj2 = subset(seurObj, subset = seurat_clusters == 0) 
#DimPlot(seurObj2, reduction = "umap")

clusters = seurObj$seurat_clusters

subsel = clusters == 0 | clusters == 1 | clusters == 4 | clusters == 5 | clusters == 9

seurObj2 = subset(seurObj, seurat_clusters == 0 | seurat_clusters == 1 | seurat_clusters == 4 | seurat_clusters == 5 | seurat_clusters == 9)
DimPlot(seurObj2, reduction = "umap")

#subsel = seurObj$seurat_clusters == 0
dsSource = seurObj$ds_source[subsel]

#now the reduced set

mergedmatsub = mergedmatFilt1[,subsel]

seurObj = CreateSeuratObject(counts = mergedmatsub, project = "merged", min.cells = 0, min.features = 0)
seurObj[["ds_source"]] <- dsSource

seurObj <- NormalizeData(seurObj, normalization.method = "LogNormalize", scale.factor = 10000)
seurObj <- FindVariableFeatures(seurObj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurObj)
seurObj <- ScaleData(seurObj, features = all.genes)
seurObj <- RunPCA(seurObj, features = VariableFeatures(object = seurObj))
seurObj <- FindNeighbors(seurObj, dims = 1:10)
seurObj <- FindClusters(seurObj, resolution = 0.5)
seurObj <- RunUMAP(seurObj, dims = 1:10)
pA = DimPlot(seurObj, reduction = "umap", group.by = "ds_source", cols=c("blue", "red"))
pA = pA + ggtitle("Batch correction - uncorrected")
pA = pA +   theme(panel.background = element_rect("white", "white", 0, 0, "white"),
                  legend.position= "bottom", legend.direction = "horizontal",#, legend.title = element_blank())
                  strip.text.x = element_text(size = 12, face = "bold"),
                  #legend.position= "none",
                  plot.title = element_text(size=14, face = "bold"),
                  strip.background = element_blank())
pA


#extract two things for calculations - the UMAP coordinates and the source
uncorrCoords = seurObj[["umap"]]@cell.embeddings
uncorrSource = dsSource
#sum(rownames(uncorrCoords) != names(uncorrSource)) #test, 0, ok


#Now, correct by down-sampling NG2 using binomial downsampling
##############################

loadStats("PBMC_NG_2")

dsc = statsPBMC_NG_2[,c(1, which(colnames(statsPBMC_NG_2) == "CPM_PBMC_NG_2_d_40"), which(colnames(statsPBMC_NG_2) == "CPM_PBMC_NG_2_d_100"))]

#so, test to use binomial downsampling
loadBug("PBMC_NG_2")
bugvNg2 = getBug("PBMC_NG_2")

cellsNotFiltered = names(dsSource)[dsSource == levels(dsSource)[2]]
bug2Filt = bugvNg2[bugvNg2$barcode %in% cellsNotFiltered,]

tmp = bug2Filt %>% group_by(gene) %>% summarise(UMIs = n(), CPM = n())
tmp$CPM = tmp$CPM*10^6/sum(tmp$CPM)


#dsManyMany = binomialDownsampling(bugvNg2, c(0.45,0.46,0.47,0.48,0.49,0.50,0.51,0.52,0.53,0.54,0.55))
vals = c(0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.40,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.50,0.51,0.52,0.53,0.54,0.55,0.9)
dsManyMany = binomialDownsampling(bug2Filt, vals)

choice_ = 25
vals[choice_] #0.45
dsMany = dsManyMany[,c(1,choice_ + 1)] #0.45 maximizes the correlation
dsManyCPM = dsMany
#mergeMany = inner_join(dsManyCPM, dsc, by="gene")
mergeMany = inner_join(dsManyCPM, tmp, by="gene")

#recalculate to CPM
for (i in 2:4) {
  mergeMany[[i]] = mergeMany[[i]] * 10^6/sum(mergeMany[[i]])
}

scale = mergeMany[[2]] / mergeMany[[4]]
NG_2mod = NG_2Mat[row.names(NG_2Mat) %in% mergeMany$gene, ]
mat = match(mergeMany$gene, row.names(NG_2mod))
sum(is.na(mat))#0, good

scaled_NG_2 = NG_2mod * scale[mat];
numNG2 = dim(scaled_NG_2)[2]

mergedmatX = merge(NGMat, scaled_NG_2, row.names(NGMat), row.names(scaled_NG_2), all.x = F, all.y = F)
mergedmatFilt2 = mergedmatX[,cellFilter]
mergedmatsub = mergedmatFilt2[,subsel]
isSecDs = dsSource == levels(dsSource)[2]

cpmNG2 = sparse_Sums(mergedmatsub[, isSecDs], T)
cpmNG = sparse_Sums(mergedmatsub[, !isSecDs], T)

cpmNG2 = cpmNG2 *10^6/sum(cpmNG2)
cpmNG = cpmNG *10^6/sum(cpmNG)

cpmNG2l = log2(cpmNG2 + 1)
cpmNGl = log2(cpmNG + 1)

getCCC(cpmNG2l, cpmNGl) # 0.9943015

#compare to original
mergedmat2 = merge(NGMat, NG_2mod, row.names(NGMat),row.names(NG_2mod), all.x = F, all.y = F)
mergedmatFilt3 = mergedmat2[,cellFilter]
mergedmatsub3 = mergedmatFilt3[,subsel]

cpmNG2__2 = sparse_Sums(mergedmatsub3[, isSecDs], T)
cpmNG__2 = sparse_Sums(mergedmatsub3[, !isSecDs], T)

cpmNG2__2 = cpmNG2__2 *10^6/sum(cpmNG2__2)
cpmNG__2 = cpmNG__2 *10^6/sum(cpmNG__2)

cpmNG2__2l = log2(cpmNG2__2 + 1)
cpmNG__2l = log2(cpmNG__2 + 1)

getCCC(cpmNG2__2l, cpmNG__2l) # 0.9913019

mergedmatFilt2 = mergedmatX[,cellFilter]
mergedmatsub = mergedmatFilt2[,subsel]

seurObj = CreateSeuratObject(counts = mergedmatsub, project = "merged", min.cells = 0, min.features = 0)
seurObj[["ds_source"]] <- dsSource
seurObj <- NormalizeData(seurObj, normalization.method = "LogNormalize", scale.factor = 10000)
seurObj <- FindVariableFeatures(seurObj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurObj)
seurObj <- ScaleData(seurObj, features = all.genes)
seurObj <- RunPCA(seurObj, features = VariableFeatures(object = seurObj))
seurObj <- FindNeighbors(seurObj, dims = 1:10)
seurObj <- FindClusters(seurObj, resolution = 0.5)
seurObj <- RunUMAP(seurObj, dims = 1:10)
pB = DimPlot(seurObj, reduction = "umap", group.by = "ds_source", cols=c("blue", "red"))
pB = pB + ggtitle("Batch correction - Corrected")
pB = pB +   theme(panel.background = element_rect("white", "white", 0, 0, "white"),
                  legend.position= "bottom", legend.direction = "horizontal",#, legend.title = element_blank())
                  strip.text.x = element_text(size = 12, face = "bold"),
                  #legend.position= "none",
                  plot.title = element_text(size=14, face = "bold"),
                  strip.background = element_blank())

pB

corrCoords = seurObj[["umap"]]@cell.embeddings
corrSource = dsSource


fig5AB = ggarrange(pA, pB, nrow=1, ncol=2, labels=c("A","B"))

ggsave(
  paste0(figure_path, "Fig5AB.png"),
  plot = fig5AB, device = "png",
  width = 9, height = 5, dpi = 300)

#######################
#Some experiments

loadStats("PBMC_NG")
loadStats("PBMC_NG_2")
s1 = getStats("PBMC_NG")
s2 = getStats("PBMC_NG_2")

s1sub = s1[,c(1, which(colnames(s1) == "FracOnes_PBMC_NG_d_100" | colnames(s1) == "UMIs_PBMC_NG_d_100" | colnames(s1) == "CPM_PBMC_NG_d_100" ))]
s1subFilt = s1sub[s1sub[[2]] >= 500, c(1,3,4)]
s2sub = s2[,c(1, which(colnames(s2) == "FracOnes_PBMC_NG_2_d_40" | colnames(s2) == "UMIs_PBMC_NG_2_d_40" | colnames(s2) == "CPM_PBMC_NG_2_d_40" ))]
s2subFilt = s2sub[s2sub[[2]] >= 500, c(1,3,4)]

jnd = inner_join(s1subFilt,s2subFilt,by="gene")
plot(jnd[[3]] ,jnd[[5]])

sort(abs(Loadings(seurObj, reduction = "pca")[, 8]), decreasing = TRUE)[1:200]

genesPC8 = names(sort(abs(Loadings(seurObj, reduction = "pca")[, 8]), decreasing = TRUE)[1:30])

genesFSCMDiff = jnd$gene[abs(jnd[[3]] - jnd[[5]]) > 0.07]
genesFSCMDiff[genesFSCMDiff %in% genesPC8]

sel2 =jnd$gene %in% genesPC8
plot(jnd[[3]][sel2] ,jnd[[5]][sel2])
lines(c(0,1),c(0,1))

sel3 = sample(dim(jnd)[1],30)
plot(jnd[[3]][sel3] ,jnd[[5]][sel3])
lines(c(0,1),c(0,1))

#test to predict the genes that differ much in FSCM
loadBug("PBMC_NG_2")
bugvNg2 = getBug("PBMC_NG_2")
cellsNotFiltered = names(dsSource)[dsSource == levels(dsSource)[2]]
bug2Filt = bugvNg2[bugvNg2$barcode %in% cellsNotFiltered,]
unPred2 = bug2Filt %>% group_by(gene) %>% summarise(UMIs = n(), CPM = n())
unPred2$CPM = unPred2$CPM*10^6/sum(unPred2$CPM)
source(paste0(sourcePath,"preseqHelpers.R"))
pred2ZTNB = upSampleAndGetMeanExprPreSeqZTNB(bug2Filt, 1000)
pred2Ds = upSampleAndGetMeanExprPreSeqDs(bug2Filt, 1000, mt=20)
pred2 = pred2Ds
predCPM2 = pred2
predCPM2[[2]] = pred2[[2]]*10^6/sum(pred2[[2]])

loadBug("PBMC_NG")
bugvNg = getBug("PBMC_NG")
cellsNotFiltered1 = names(dsSource)[dsSource == levels(dsSource)[1]]
bugFilt = bugvNg[bugvNg$barcode %in% cellsNotFiltered1,]
unPred = bugFilt %>% group_by(gene) %>% summarise(UMIs = n(), CPM = n())
unPred$CPM = unPred$CPM*10^6/sum(unPred$CPM)
predZTNB = upSampleAndGetMeanExprPreSeqZTNB(bugFilt, 1000/0.45)
predDs = upSampleAndGetMeanExprPreSeqDs(bugFilt, 1000/0.45, mt=20)
pred = predDs
predCPM = pred
predCPM[[2]] = pred[[2]]*10^6/sum(pred[[2]])

unPredicted = inner_join(unPred,unPred2, by="gene")
unPredicted[unPredicted$gene %in%  genesPC8,]
predicted = inner_join(predCPM,predCPM2, by="gene")
predicted[predicted$gene %in%  genesPC8,]

genesPC8Inter = genesFSCMDiff[genesFSCMDiff %in% genesPC8]

unPredFilt = unPredicted[unPredicted$gene %in% genesPC8Inter,]
unPredFilt
predFilt = predicted[predicted$gene %in% genesPC8Inter,]
predFilt

plotHists = function(b1,b2,gene) {
  gbug1 = b1[b1$gene == gene,]
  gbug2 = b2[b2$gene == gene,]
  hist(gbug1$count, breaks=seq(0.5, max(gbug1$count)+0.5, by=1), main = paste0("b1: ",gene))
  hist(gbug2$count, breaks=seq(0.5, max(gbug2$count)+0.5, by=1), main = paste0("b2: ",gene))
}
plotHists(bugvNg, bugvNg2, "AC010970.1")
plotHists(bugvNg, bugvNg2, "MTATP6P1")
plotHists(bugvNg, bugvNg2, "RPL26")
plotHists(bugvNg, bugvNg2, "RPL7")
plotHists(bugvNg, bugvNg2, "RPS10")
plotHists(bugvNg, bugvNg2, "RPS2")
gbug1 = bugvNg[bugvNg$gene == "AC010970.1",]
gbug2 = bugvNg2[bugvNg2$gene == "AC010970.1",]
hist(gbug1$count, breaks=seq(0.5, max(gbug1$count)+0.5, by=1))
hist(gbug2$count, breaks=seq(0.5, max(gbug2$count)+0.5, by=1))


#Calculate the 10 nearest neighbors for each cell


getKnnFrac = function(crds, src, numNeighbors) {
  nn = length(src)
  xs = as.data.frame(crds)$UMAP_1
  ys = as.data.frame(crds)$UMAP_2
  res = rep(NA, nn)
  for (i in 1:nn) {
    sci = src[i]
    sc = src[-i]
    dists = sqrt((xs[i] - xs[-i])^2 + (ys[i] - ys[-i])^2) #remove self
    srt = sort(dists, index.return=TRUE)
    res[i] = sum(sci == sc[srt$ix[1:numNeighbors]])/numNeighbors
  }
  return (mean(res))
}

getKnnFrac(uncorrCoords, uncorrSource, 10)#0.8185659
getKnnFrac(corrCoords, corrSource, 10)#0.6765698

#test case:
#     0 1 2.2
#0    0 0 0
#0.9  1 1 0
#2    1 1 0
#3.1  1 1 0
testCoords = tibble(UMAP_1 = rep(c(0,1,2.2),4), UMAP_2 = c(rep(0,3), rep(0.9,3), rep(2,3), rep(3.1,3)))
testSrc = c(0,0,0,1,1,0,1,1,0,1,1,0)
expRes = mean(c(0.5, 0.5, 1, 0.5, 0.5, 1, 1, 1, 1, 1, 1, 0.5))
res = getKnnFrac(testCoords, testSrc, 2)
expRes - res < 0.0000001 #ok



################
# test the binomial downsampling function
#################

#loadBug("PBMC_V3_3")
#bugv3_3 = getBug("PBMC_V3_3")
#bug = bugv3_3
#t = 0.1
#ds = binomialDownsampling(bug, c(0.1,0.5))


