sourcePath = "C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/"

source(paste0(sourcePath,"ButterflyHelpers.R"))

localDataPath = "C:/Work/R/ButterflyQuant/"
cachedDataPath = paste0(localDataPath,"savedData/")
dataPath = "E:/Butterfly/"
figure_data_path = "C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/FigureData/"
figure_path = "Z:/projects/Butterfly/figures/"

library(ggplot2)
library("ggpubr")
library("tidyverse")
library("dplyr")

loadStats("EVALPBMC_DS")
loadStats("EVALPBMC")
loadStats("PBMC_V3")
loadStats("PBMC_V3_2")
loadStats("PBMC_V3_3")
loadStats("PBMC_V2")
loadStats("PBMC_NG")
loadStats("PBMC_NG_2")

umisList = list(tibble(gene=statsEVALPBMC$gene, umis=statsEVALPBMC$UMIs_EVALPBMC_d_100), 
                tibble(gene=statsPBMC_V3$gene, umis=statsPBMC_V3$UMIs_PBMC_V3_d_100), 
                tibble(gene=statsPBMC_V3_2$gene, umis=statsPBMC_V3_2$UMIs_PBMC_V3_2_d_100), 
                tibble(gene=statsPBMC_V3_3$gene, umis=statsPBMC_V3_3$UMIs_PBMC_V3_3_d_100), 
                tibble(gene=statsPBMC_V2$gene, umis=statsPBMC_V2$UMIs_PBMC_V2_d_100), 
                tibble(gene=statsPBMC_NG$gene, umis=statsPBMC_NG$UMIs_PBMC_NG_d_100), 
                tibble(gene=statsPBMC_NG_2$gene, umis=statsPBMC_NG_2$UMIs_PBMC_NG_2_d_100))

fracOnesList = list(tibble(gene = statsEVALPBMC$gene, fracOnes = statsEVALPBMC$FracOnes_EVALPBMC_d_100),
                    tibble(gene = statsPBMC_V3$gene, fracOnes = statsPBMC_V3$FracOnes_PBMC_V3_d_100),
                    tibble(gene = statsPBMC_V3_2$gene, fracOnes = statsPBMC_V3_2$FracOnes_PBMC_V3_2_d_100),
                    tibble(gene = statsPBMC_V3_3$gene, fracOnes = statsPBMC_V3_3$FracOnes_PBMC_V3_3_d_100),
                    tibble(gene = statsPBMC_V2$gene, fracOnes = statsPBMC_V2$FracOnes_PBMC_V2_d_100),
                    tibble(gene = statsPBMC_NG$gene, fracOnes = statsPBMC_NG$FracOnes_PBMC_NG_d_100),
                    tibble(gene = statsPBMC_NG_2$gene, fracOnes = statsPBMC_NG_2$FracOnes_PBMC_NG_2_d_100))


umis = umisList[[1]]
fos = fracOnesList[[1]]
for (i in 2:length(umisList)) {
  umis = inner_join(umis, umisList[[i]], by="gene")
  fos = inner_join(fos, fracOnesList[[i]], by="gene")
}

mns = rep(0,dim(fos)[1])

#calculate weighted mean with a loop
for (i in 1:(dim(fos)[1])) {
  mns[i] = weighted.mean(fos[i,-1], umis[i,-1])
}

library(matrixStats)
mns = rowWeightedMeans(fos[,-1], w = NULL)

#test: the numbers look reasonably right, didn't check any details
#so, the weighted mean gives a reasonably good estimate of which genes that are most extreme
#there are of course error factors such as different saturation levels for datasets, but we just ignore that here

srt = sort(mns, index.return=T, decreasing = T)

srtumis = umis[srt$ix,]
srtumis

totumis = rowSums(srtumis[,-1])
totumis #looks like there is no need for filtering due to lack of data, they all have several hundred umis

srtumis[1:100,]
srt$x[1:100]

srtumis$gene[1:100]

#now write the genes to file
write_csv2(srtumis[1:200,1], paste0(figure_data_path, "extreme_genes.txt"))

#get the number for NEUROD1, also in drop-seq
paste0("Fraction of single-copy molecules in 10x data for NEUROD1: ", mns[umis$gene == "NEUROD1"])#0.897550111358575
statsEVALPBMC_DS$UMIs_EVALPBMC_DS_d_100[statsEVALPBMC_DS$gene == "NEUROD1"]
dsFos = statsEVALPBMC_DS$FracOnes_EVALPBMC_DS_d_100[statsEVALPBMC_DS$gene == "NEUROD1"]
paste0("Fraction of single-copy molecules in Drop-Seq data for NEUROD1: ", dsFos)#0.819548872180451

#some experimental code below
#see if we can find an interesting gene which differs much between dropseq and 10x:
#genes = srtumis$gene[1:200]
#statsEVALPBMC_DSFilt = statsEVALPBMC_DS[statsEVALPBMC_DS$UMIs_EVALPBMC_DS_d_100 > 100,]
  
#statsEVALPBMC_DSFiltSel = statsEVALPBMC_DSFilt[statsEVALPBMC_DSFilt$gene %in% genes,]
#statsEVALPBMC_DSFiltSel2 = statsEVALPBMC_DSFiltSel[statsEVALPBMC_DSFiltSel$FracOnes_EVALPBMC_DS_d_100 < 0.70, ]
#this is not very useful...

#ACE2, the corona gene, is not in those datasets (perhaps lost in inner_join). It should be present in the lung cancer though:

#loadStats("LC")
#loadStats("EVALPBMC_SW")
#loadStats("EVALPBMC_DS")
#loadStats("LC")

#statsLC[statsLC$gene == "ACE2",]
#statsEVALPBMC[statsEVALPBMC$gene == "ACE2",]
#statsEVALPBMC_SW[statsEVALPBMC_SW$gene == "ACE2",]
#statsEVALPBMC_DS[statsEVALPBMC_DS$gene == "ACE2",]

#statsPBMC_V3[statsPBMC_V3$gene == "ACE2",]
#statsPBMC_V3_2[statsPBMC_V3_2$gene == "ACE2",]
#statsPBMC_V3_3[statsPBMC_V3_3$gene == "ACE2",]
#statsPBMC_V2[statsPBMC_V2$gene == "ACE2",]
#statsPBMC_NG[statsPBMC_NG$gene == "ACE2",]
#statsPBMC_NG_2[statsPBMC_NG_2$gene == "ACE2",]


