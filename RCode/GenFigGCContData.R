#
# Generates the data for the plots in Fig. 1
#

#before anything else, you need to setup paths:
#for example
#source("C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/paths.R")
#or, in colab:
#source("GRNP_2020/RCode/pathsGoogleColab.R")

source(paste0(sourcePath,"ButterflyHelpers.R"))
#source(paste0(sourcePath,"PreseqHelpers.R"))

library(tidyverse)

###############
# Gene length
##############
#BiocManager::install("GenomicFeatures")
library("GenomicFeatures")

txdbM = makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", transcript_ids=NULL, circ_seqs=DEFAULT_CIRC_SEQS, filter=NULL, id_prefix="ensembl_", host="www.ensembl.org", port=80, taxonomyId=NA, miRBaseBuild=NA)

tlM = transcriptLengths(txdbM, with.cds_len=FALSE)
colnames(tlM)[c(2,3)] = c("tx","gene") #rename to match with the other datasets



library(biomaRt)


#get versions
listMarts()#version 98
ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(ensembl)

##############
#GC Content
##############

#load genome for mus musculus
#BiocManager::install("BSgenome")
#BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
library(BSgenome.Mmusculus.UCSC.mm10)
genomeM <- BSgenome.Mmusculus.UCSC.mm10
#BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdbM <- TxDb.Mmusculus.UCSC.mm10.knownGene
transcriptsM <- exonsBy(txdbM, by="tx", use.names=TRUE)


seqsM = extractTranscriptSeqs(genomeM, transcriptsM)

#calculate gc content in the strings

gcContent <- function(x) {
  return (letterFrequency(x, c("GC"), OR="|", as.prob=TRUE))
} 


#test TC002 - gcContent
b = BString("GGCCGA")
gcContent(b)#should be 5/6 = 0.8333333, ok!


gcFullLength = gcContent(seqsM)
txs = seqsM@ranges@NAMES
#need to remove the version from the transcript name
txs = substr(txs, 1, 18)

#these have transcript ids only, not gene ids. So, merge with the length to get gene id
gcs = tibble(tx = txs, gc = gcFullLength[,1])

#convert the genes
tr2g = read.table(paste0(dataPath,"EVAL/bus_output/transcripts_to_genes.txt"), stringsAsFactors = F)
lookupTable = tr2g[,2:3]
lookupTable= unique(lookupTable)
lookupTable[[1]] = substr(lookupTable[[1]], 1, str_length("ENSMUSG00000087582"))

library(qdapTools)
outGenes = lookup(tlM$gene, lookupTable)
#length(outGenes)#140725
#dim(tlM) #140725      5, ok
tlM$gene = outGenes

gcsMMerged = inner_join(gcs, tlM, by="tx")
#remove NAs (genes that don't have a gene name)
gcsMMergedFilt = gcsMMerged[!is.na(gcsMMerged$gene),]

#take mean of all transcripts for each gene
gctl = gcsMMergedFilt %>% group_by(gene) %>% summarize(gc=mean(gc), txlen = mean(tx_len))

#test
#gcsMMergedFilt[gcsMMergedFilt$gene == "Pbsn",]
#1 ENSMUST00000000003 0.401 137139 Pbsn      7    902
#2 ENSMUST00000114041 0.396 137140 Pbsn      6    697
#txlen should be (902+697)/2 = 799.5
#gc should be (0.401 + 0.396)/2 = 0.3985
#real values:
#gctl[gctl$gene == "Pbsn", ]
#Pbsn 0.399  800., ok (rounded)

#Get the data - we get FSCM and CU from the full EVAL dataset:
loadStats("EVAL")
stats = getStats("EVAL")
d = stats[,c(1,which((colnames(stats) == "FracOnes_EVAL_d_100") | (colnames(stats) == "CountsPerUMI_EVAL_d_100") | (colnames(stats) == "UMIs_EVAL_d_100")))]
colnames(d)[2:4] = c("UMIs","FSCM", "CU")
d = d[d$UMIs >= 30,]

#merge everything together
gctlcu = inner_join(gctl, d, by="gene")
dim(gctlcu)#11826

saveRDS(gctlcu, paste0(figure_data_path, "gc.RDS"))
