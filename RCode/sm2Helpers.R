library(tidyverse)
#library(BUSpaRse)
library(textTinyR)
library(DescTools)
library(qdapTools)


ReadKallistoSm2TxMatrix = function(tr2gFile, trLengthsFile, outputDir, fnMatrix="matrix.abundance.mtx", fnTx="transcripts.txt", fnCells="matrix.cells") {
  tr2g = read_tsv(tr2gFile, col_names = FALSE)
  colnames(tr2g) = c("Tx", "ensembl", "symbol")
  
  #so, the first line is the number of rows and columns
  dsMatTmp = read_tsv(paste0(outputDir, fnMatrix), col_names = FALSE, comment="%%")
  firstLine = dsMatTmp[1,]
  dsMat = dsMatTmp[-1,]
  ds = sparseMatrix(i= dsMat[[2]], j = dsMat[[1]], x = dsMat[[3]], dims = c(firstLine[[2]], firstLine[[1]]))
  #unsparse for practical reasons
  ds = as.matrix(ds)

  cellIds = read_tsv(paste0(outputDir, fnCells), col_names = FALSE)
  txs = read_tsv(paste0(outputDir, fnTx), col_names = FALSE)
  rownames(ds) = txs[[1]]
  colnames(ds) = cellIds[[1]]
  
  dsTPM = ds

  #get transcript lengths, join the datasets, divide by tx length and scale to 10^6
  dsLengths = read_tsv(trLengthsFile)[,c(1,2)] #use effective length
  #colnames(dsLengths) = c("Tx", "Length")
  colnames(dsLengths) = c("x", "y")
  
  lengths = lookup(rownames(dsTPM), dsLengths)
  
  
  #loop through the cells to be sure they are handled correctly
  for (i in 1:ncol(dsTPM)) {
    dsTPM[,i] = dsTPM[,i]/lengths
    dsTPM[,i] = dsTPM[,i] * 10^6/sum(dsTPM[,i])
  }
  
  return (list(ds, dsTPM))
}
