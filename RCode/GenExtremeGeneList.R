#before anything else, you need to setup paths:
#for example
#source("C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/paths.R")
#or, in colab:
#source("GRNP_2020/RCode/pathsGoogleColab.R")

source(paste0(sourcePath,"ButterflyHelpers.R"))

library(readr)

#######################
# Human 10X
#######################

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

#test: the numbers look reasonably right
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

toWrite = tibble(gene=srtumis$gene, UMIs = totumis, FSCM=srt$x)

#now write the genes to file
write_delim(toWrite, paste0(figure_path, "human_10x_FSCM_genes.txt"), "\t")

#get the number for NEUROD1, also in drop-seq
paste0("Fraction of single-copy molecules in 10x data for NEUROD1: ", mns[umis$gene == "NEUROD1"])#0.897550111358575
statsEVALPBMC_DS$UMIs_EVALPBMC_DS_d_100[statsEVALPBMC_DS$gene == "NEUROD1"]
dsFos = statsEVALPBMC_DS$FracOnes_EVALPBMC_DS_d_100[statsEVALPBMC_DS$gene == "NEUROD1"]
paste0("Fraction of single-copy molecules in Drop-Seq data for NEUROD1: ", dsFos)#0.819548872180451

#########################
# Human Drop-Seq
#########################
loadStats("EVALPBMC_DS")

toWrite = tibble(gene=statsEVALPBMC_DS$gene, UMIs = statsEVALPBMC_DS$UMIs_EVALPBMC_DS_d_100, FSCM=statsEVALPBMC_DS$FracOnes_EVALPBMC_DS_d_100)
#filter on 20, otherwise it is pretty random
toWrite = toWrite[toWrite$UMIs >= 20, ]
srt = sort(toWrite$FSCM, index.return=TRUE, decreasing=TRUE)
toWrite = toWrite[srt$ix,]

#now write the genes to file
write_delim(toWrite, paste0(figure_path, "human_Drop-Seq_FSCM_genes.txt"), "\t")

#######################
# Mouse 10X
#######################

loadStats("EVAL")
loadStats("MRET2")

umisList = list(tibble(gene=statsEVAL$gene, umis=statsEVAL$UMIs_EVAL_d_100), 
                tibble(gene=statsMRET2$gene, umis=statsMRET2$UMIs_MRET2_d_100))

fracOnesList = list(tibble(gene = statsEVAL$gene, fracOnes = statsEVAL$FracOnes_EVAL_d_100),
                    tibble(gene = statsMRET2$gene, fracOnes = statsMRET2$FracOnes_MRET2_d_100))


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

srt = sort(mns, index.return=T, decreasing = T)

srtumis = umis[srt$ix,]
srtumis

totumis = rowSums(srtumis[,-1])
totumis #looks like there is no need for filtering due to lack of data, they all have several hundred umis

toWrite = tibble(gene=srtumis$gene, UMIs = totumis, FSCM=srt$x)

#now write the genes to file
write_delim(toWrite, paste0(figure_path, "mouse_10x_FSCM_genes.txt"), "\t")

#########################
# Mouse Drop-Seq
#########################
loadStats("MRET")

toWrite = tibble(gene=statsMRET$gene, UMIs = statsMRET$UMIs_MRET_d_100, FSCM=statsMRET$FracOnes_MRET_d_100)
#filter on 20, otherwise it is pretty random
toWrite = toWrite[toWrite$UMIs >= 20, ]
srt = sort(toWrite$FSCM, index.return=TRUE, decreasing=TRUE)
toWrite = toWrite[srt$ix,]

#now write the genes to file
write_delim(toWrite, paste0(figure_path, "mouse_Drop-Seq_FSCM_genes.txt"), "\t")






