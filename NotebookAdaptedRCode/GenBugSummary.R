sourcePath = "GRNP_2020/NotebookAdaptedRCode/"

source(paste0(sourcePath,"ButterflyHelpers.R"))


#used for investigating single-copy molecules
ClosestDists <- function(bug, subBug, UMILength) {
  #make sure subbug is max 2000 UMIs to shorten the execution time:
  if (dim(subBug)[1] > 2000) {
    print("Down-sampling to 2000 UMIs")
    samp = sort(sample(1:(dim(subBug)[1]), 2000))
    subBug = subBug[samp,]
  }

  currCell = ""
  dists = rep(0,UMILength)
  print(paste0("Will process ", dim(subBug)[1], " UMIs"))
  for (i in 1:(dim(subBug)[1])) {
    if (i%%100 == 0) {
      print(i)
    }
    if (subBug$barcode[i] != currCell) {
      currCell = subBug$barcode[i];
      currCellBug = bug[bug$barcode==currCell,]
    }
    currUMI = subBug$UMI[[i]]
    bestDist = UMILength;
    for (j in 1:(dim(currCellBug)[1])) {
      dist = stringdist::stringdist(currCellBug$UMI[[j]], currUMI, method = "hamming")
      if (dist > 0 & dist < bestDist) {
        bestDist = dist;
        if (bestDist == 1){
          break
        }
      }
    }
    dists[bestDist] = dists[bestDist] + 1;
  }
  return (dists)
}


#highSingleCopyUMIGene can be Vmn1r13 for mouse
#lowSingleCopyUMIGene can be Ubb for mouse
genBugSummary <- function(dsid, lowSingleCopyUMIGene, highSingleCopyUMIGene, UMILen, fig_data_path = figure_data_path) {
  loadBug(dsid, fig_data_path = fig_data_path)
  bug = get(paste0("bug",dsid), envir=.GlobalEnv)
  
  #Some simple statistics
  totCells = length(unique(bug$barcode))
  totUMIs = dim(bug)[1]
  totCounts = sum(bug$count)
  countsPerUMI = mean(bug$count)
  h1 = hist(bug$count, breaks=seq(0.5, max(bug$count)+0.5, by=1), plot=F)
  UMIsPerCell = dim(bug)[1]/totCells
  countsPerCell = totCounts/length(unique(bug$barcode))
  totFracOnes = h1$density[1]

  #compare the fraction of close barcodes a gene with many single-copy molecules compared to one with few
  #the purpose is that most false molecules from read errors etc. will mostly manifest as single-copy molecules
  #Vmn1r13, almost only ones
  genebug = bug[bug$gene==highSingleCopyUMIGene,]
  hscgPosCells = length(unique(genebug$barcode))
  hscgUMIs = length(genebug$barcode)
  distsHscg = ClosestDists(bug, genebug, UMILen)
  # investigate if the UMI barcodes are closer to one another for genes 
  # with a higher single-copy molecule fraction
  genebug = bug[bug$gene==lowSingleCopyUMIGene,]
  lscgPosCells = length(unique(genebug$barcode))
  lscgUMIs = length(genebug$barcode)
  distsLscg = ClosestDists(bug, genebug, UMILen)
  
  distsHscgStrings = paste(as.character(distsHscg), collapse=", ")
  distsLscgStrings = paste(as.character(distsLscg), collapse=", ")
  
  distsHscgFrac = distsHscg/sum(distsHscg)
  distsHscgFracStrings = paste(as.character(distsHscgFrac), collapse=", ")
  distsLscgFrac = distsLscg/sum(distsLscg)
  distsLscgFracStrings = paste(as.character(distsLscgFrac), collapse=", ")
  
  # investigate if the UMI barcodes from single-copy molecules 
  # have closer neighbors than for molecules with more copies
  bugsel = bug

  bug1cpy = bugsel[bugsel$count==1,]
  dists1cpy = ClosestDists(bugsel, bug1cpy, UMILen)
  
  bug2cpy = bugsel[bugsel$count==2,]
  dists2cpy = ClosestDists(bugsel, bug2cpy, UMILen)
  
  bug3cpy = bugsel[bugsel$count>=3,]
  dists3cpy = ClosestDists(bugsel, bug3cpy, UMILen)

  dists1cpyStrings = paste(as.character(dists1cpy), collapse=", ")
  dists1cpyFrac = dists1cpy/sum(dists1cpy)
  dists1cpyFracStrings = paste(as.character(dists1cpyFrac), collapse=", ")
  
  dists2cpyStrings = paste(as.character(dists2cpy), collapse=", ")
  dists2cpyFrac = dists2cpy/sum(dists2cpy)
  dists2cpyFracStrings = paste(as.character(dists2cpyFrac), collapse=", ")
  
  dists3cpyStrings = paste(as.character(dists3cpy), collapse=", ")
  dists3cpyFrac = dists3cpy/sum(dists3cpy)
  dists3cpyFracStrings = paste(as.character(dists3cpyFrac), collapse=", ")
  
  #now write
  fileConn<-file(paste0(fig_data_path, dsid, "/ds_summary.txt"))
  writeLines(c(paste0("Dataset: ", dsid, "\n"),
             paste0("totUMIs: ", totUMIs),
             paste0("totCells: ", totCells),
             paste0("totCounts: ", totCounts),
             paste0("countsPerUMI: ", countsPerUMI),
             paste0("UMIsPerCell: ", UMIsPerCell),
             paste0("countsPerCell: ", countsPerCell),
             paste0("totFracOnes: ", totFracOnes),
             paste0("FracMolWithUMIDistToNeighborH: ", distsHscgStrings),
             paste0("FracMolWithUMIDistToNeighborL: ", distsLscgStrings),
             paste0("FracMolWithUMIDistToNeighborHFrac: ", distsHscgFracStrings),
             paste0("FracMolWithUMIDistToNeighborLFrac: ", distsLscgFracStrings),
             paste0("FracMolWithUMIDistToNeighbor1cpy: ", dists1cpyStrings),
             paste0("FracMolWithUMIDistToNeighbor2cpy: ", dists2cpyStrings),
             paste0("FracMolWithUMIDistToNeighbor>=3cpy: ", dists3cpyStrings),
             paste0("FracMolWithUMIDistToNeighbor1cpyFrac: ", dists1cpyFracStrings),
             paste0("FracMolWithUMIDistToNeighbor2cpyFrac: ", dists2cpyFracStrings),
             paste0("FracMolWithUMIDistToNeighbor>=3cpyFrac: ", dists3cpyFracStrings)
  ), fileConn)
  close(fileConn)
}

