#before anything else, you need to setup paths:
#for example
#source("C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/paths.R")
#or, in colab:
#source("GRNP_2020/RCode/pathsGoogleColab.R")

source(paste0(sourcePath, "ButterflyHelpers.R"))


#library(stringr)
library(dplyr)
library(qdapTools)

#convert gene indices to gene symbols
geneIndices2Symbols <- function(geneIndices, genesFile, Tr2gFile) {
  genenames = read.table(genesFile, stringsAsFactors = F)$V1
  tr2g = read.table(Tr2gFile, stringsAsFactors = F)
  lookupTable = tr2g[,2:3]
  lookupTable= unique(lookupTable)
#  lookupTableTranscr = tr2g[,c(1,3)]
  numGenes = length(geneIndices)

  inGenes = geneIndices;#allocate

  for (i in 1:numGenes) {
    inGenes[i] = genenames[1 + as.numeric(geneIndices[i])]
  }
  #check that the ens genes are really unique in the lookup table:
#  length(lookupTable$V2)
#  length(unique(lookupTable$V2)) #ok

  outGenes = lookup(inGenes, lookupTable)

  return (outGenes)
}

#dir should include slash at the end
#conserveMem will make the removal of multimapped reads approx. 4 times slower, but will preserve memory, which is good for notebooks
#when running large datasets
readBug <- function(dir, conserveMem = F) {
  print(paste0("Reading BUG from ", dir, " ..."))
  bug = read.table(paste0(dir,"bus_output/bug.txt"), stringsAsFactors = F)
  print("Filtering multi-mapped reads...")
  
  if (!conserveMem) {
    isMult = sapply(bug[,3],function(s) grepl(",",s, fixed=T))
  } else {
    #do this with a for loop instead - the extra execution time is negligible, and we get more control of the memory use (notebooks run out of memory here for the LC dataset)
    print("Using a slower processing method to preserve memory...")
    sz = dim(bug)[1]
    isMult = logical(sz)
    for (i in 1:sz) {
  	  isMult[i] = grepl(",",bug[i,3], fixed=T) 
    }
  }

  print (paste0("Fraction multi-mapped reads: ", sum(isMult) / dim(bug)[1]))

  uniquelymapped = bug[!isMult,] #just ignore the multimapped reads
  rm(bug)

  print("Converting genes...")

  #now convert gene index to gene id
  uniquelymapped$V3 = geneIndices2Symbols(uniquelymapped$V3, paste0(dir,"bus_output/coll.genes.txt"), paste0(dir,"bus_output/transcripts_to_genes.txt"))

  colnames(uniquelymapped) = c("barcode", "UMI", "gene", "count")

  print("Done")

  return (uniquelymapped)
}


formatFracName<- function(name, frac, type) {
  paste0(type, "_", name,"_d_", frac*100)
}

fracOnesFunc <- function(d) {
  sum(d == 1) / length(d)
}


createStandardBugsData <- function(bugdir, name, fracs, UmisPerCellLimit = 200, fig_data_path = figure_data_path, conserveMem = F) {
  #Generate data
  print(paste0("Generating data for ", name))
  bug = readBug(bugdir, conserveMem)
  #filter out low quality cells
  #Should have more than 200 UMIs
  UMIsPerCell = bug %>% group_by(barcode) %>% tally()
  sum(UMIsPerCell$n > UmisPerCellLimit)
  filtBug = bug[bug$barcode %in% UMIsPerCell$barcode[UMIsPerCell$n > UmisPerCellLimit],]
  rm(bug)
  
  #skip mitochondrial content for now; don't think it matters for this application
  #should have mitochondrial content of less than 10%
  #genelist = filtBug %>% group_by(barcode) %>% do(genes=c(.$gene))
  
  #mc = rep(0,length(genelist$barcode))
  
  #for (i in 1:length(genelist)) {
  #  tmp <- grep(pattern = "^mt-", x = genelist$genes[[i]], value = TRUE)
  #  mc[i] = length(tmp) /length(genelist$genes[[i]])
  #}
  
  #create folder
  subpath = paste0(fig_data_path, name, "/")
  if (!file.exists(paste0(fig_data_path, name))){ #there must be no slash at the end here...
    dir.create(file.path(fig_data_path, name))
  }
  #generate proper variable names before saving
  fracsName = paste0("fracs",name)
  statsName = paste0("stats",name)
  
  
  #downsample bug, save it and collect stats
  statsList = vector(mode = "list", length = length(fracs))
  print(paste0("Down-sampling in total ", length(fracs), " bugs:"))
  for (i in 1:length(fracs)) {
    print(paste0(i, ": Down-sampling to ", fracs[[i]]))
    if (fracs[[i]] == 1) {
      dsBug = filtBug
    } else {
      dsBug = downSampleBUG(filtBug,fracs[[i]])
    }
    #save the bug
    print("saving BUG...")
    bugName = genBugObjName(name, fracs[[i]])
    bugFileName = genBugFileName(fracs[[i]])
    namedSave(list(dsBug), list(bugName), paste0(subpath, bugFileName))
    
    #create stats
    tmp = dsBug %>% group_by(gene) %>% summarise(UMIs = n(), counts=sum(count), CPM = n(), fracOnes = fracOnesFunc(count), countsPerUMI = mean(count))
    #CPM needs to be fixed, only UMI counts right now
    tmp$CPM = tmp$CPM*10^6/sum(tmp$CPM)
    #set the right column names here
    cn = c("gene",
           formatFracName(name, fracs[[i]], "UMIs"),
           formatFracName(name, fracs[[i]], "Counts"),
           formatFracName(name, fracs[[i]], "CPM"),
           formatFracName(name, fracs[[i]], "FracOnes"),
           formatFracName(name, fracs[[i]], "CountsPerUMI")
    )
    colnames(tmp) = cn
    statsList[[i]] = tmp
	
	rm(dsBug)
    
  }
  print("Done")
  
  #now merge stats to one table
  stats = statsList[[1]]
  if (length(fracs) > 1) {
    for(i in 2:length(fracs)) {
      stats = inner_join(stats, statsList[[i]], by="gene")
    }
  }
  
  print("Saving stats...")
  namedSave(list(stats), list(statsName), paste0(subpath, "Stats.RData"))
}


