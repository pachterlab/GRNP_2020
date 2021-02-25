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

ggsave(
  paste0(figure_path, "FigS_batch.png"),
  plot = supFigA, device = "png",
  width = 8, height = 5, dpi = 300)



#seurObj2 = subset(seurObj, subset = seurat_clusters == 0 | seurat_clusters == 3 | seurat_clusters == 5 | seurat_clusters == 6) 
#seurObj2 = subset(seurObj, subset = seurat_clusters == 0) 
#DimPlot(seurObj2, reduction = "umap")

clusters = seurObj$seurat_clusters

subsel = clusters == 0 | clusters == 3 | clusters == 5 | clusters == 6
#subsel = seurObj$seurat_clusters == 0
dsSource = seurObj$ds_source[subsel]

#now the reduced set

mergedmatsub = mergedmat[,subsel]

seurObj = CreateSeuratObject(counts = mergedmatsub, project = "merged", min.cells = 0, min.features = 0)
seurObj[["percent.mt"]] <- PercentageFeatureSet(seurObj, pattern = "^MT-")
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


#Now, correct by down-sampling NG2 to 0.49
##############################

loadStats("PBMC_NG_2")

dsc = statsPBMC_NG_2[,c(1, which(colnames(statsPBMC_NG_2) == "CPM_PBMC_NG_2_d_40"), which(colnames(statsPBMC_NG_2) == "CPM_PBMC_NG_2_d_100"))]

#so, test to use binomial downsampling
loadBug("PBMC_NG_2")
bugvNg2 = getBug("PBMC_NG_2")

dsManyMany = binomialDownsampling(bugvNg2, c(0.45,0.46,0.47,0.48,0.49,0.50,0.51,0.52,0.53,0.54,0.55))

dsMany = dsManyMany[,c(1,6)] #0.49 maximizes the correlation
dsManyCPM = dsMany
mergeMany = inner_join(dsManyCPM, dsc, by="gene")

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

cpmNG2 = sparse_Sums(mergedmatX[, 1:numNG2], T)
cpmNG = sparse_Sums(mergedmatX[, (numNG2+1):(dim(mergedmatX)[2])], T)

cpmNG2 = cpmNG2 *10^6/sum(cpmNG2)
cpmNG = cpmNG *10^6/sum(cpmNG)

cpmNG2l = log2(cpmNG2 + 1)
cpmNGl = log2(cpmNG + 1)

getCCC(cpmNG2l, cpmNGl) # 0.9935431

#compare to original
mergedmat2 = merge(NG_2mod, NGMat,  row.names(NG_2mod), row.names(NGMat),all.x = F, all.y = F)

cpmNG2__2 = sparse_Sums(mergedmat2[, 1:numNG2], T)
cpmNG__2 = sparse_Sums(mergedmat2[, (numNG2+1):(dim(mergedmat2)[2])], T)

cpmNG2__2 = cpmNG2__2 *10^6/sum(cpmNG2__2)
cpmNG__2 = cpmNG__2 *10^6/sum(cpmNG__2)

cpmNG2__2l = log2(cpmNG2__2 + 1)
cpmNG__2l = log2(cpmNG__2 + 1)

getCCC(cpmNG2__2l, cpmNG__2l) # 0.9904393

mergedmatsub = mergedmatX[,subsel]

seurObj = CreateSeuratObject(counts = mergedmatsub, project = "merged", min.cells = 0, min.features = 0)
#seurObj[["percent.mt"]] <- PercentageFeatureSet(seurObj, pattern = "^MT-")
seurObj[["ds_source"]] <- dsSource
seurObj <- NormalizeData(seurObj, normalization.method = "LogNormalize", scale.factor = 10000)
seurObj <- FindVariableFeatures(seurObj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurObj)
seurObj <- ScaleData(seurObj, features = all.genes)
#seurObj = SCTransform(seurObj, vars.to.regress = "percent.mt", verbose = TRUE)
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

