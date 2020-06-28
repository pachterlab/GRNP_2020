sourcePath = "GRNP_2020/NotebookAdaptedRCode/"
source(paste0(sourcePath, "paths.R"))
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
readBug <- function(dir) {
  print(paste0("Reading BUG from ", dir, " ..."))
  bug = read.table(paste0(dir,"bus_output/bug.txt"), stringsAsFactors = F)
  print("Filtering multi-mapped reads...")
  isMult = sapply(bug[,3],function(s) grepl(",",s, fixed=T))

  print (paste0("Fraction multi-mapped reads: ", sum(isMult) / dim(bug)[1]))

  uniquelymapped = bug[!isMult,] #just ignore the multimapped reads

  print("Converting genes...")

  #now convert gene index to gene id
  uniquelymapped$V3 = geneIndices2Symbols(uniquelymapped$V3, paste0(dir,"bus_output/coll.genes.txt"), paste0(dir,"bus_output/transcripts_to_genes.txt"))

  colnames(uniquelymapped) = c("barcode", "UMI", "gene", "count")

  print("Done")

  return (uniquelymapped)
}

downsampleBugs <- function(bug, fracs) {
  print(paste0("Down-sampling in total ", length(fracs), " bugs:"))
  bugs = vector(mode = "list", length = length(fracs))
  for (i in 1:length(fracs)) {
    print(paste0(i, ": Down-sampling to ", fracs[[i]]))
    if (fracs[[i]] == 1) {
      bugs[[i]] = bug
    } else {
      bugs[[i]] = downSampleBUG(bug,fracs[[i]])
    }
  }
  print("Done")
  return(bugs)
}

formatFracName<- function(name, frac, type) {
  paste0(type, "_", name,"_d_", frac*100)
}

fracOnesFunc <- function(d) {
  sum(d == 1) / length(d)
}


getStatsFromBugs <- function(dsBugs, fracs, name) {
  stats = vector(mode = "list", length = length(fracs))
  for (i in 1:length(dsBugs)) {
    tmp = dsBugs[[i]] %>% group_by(gene) %>% summarise(UMIs = n(), counts=sum(count), CPM = n(), fracOnes = fracOnesFunc(count), countsPerUMI = mean(count))
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
    stats[[i]] = tmp
  }

  #now merge to one table
  stat = stats[[1]]
  if (length(fracs) > 1) {
    for(i in 2:length(fracs)) {
      stat = inner_join(stat, stats[[i]], by="gene")
    }
  }

  return (stat)
}

createStandardBugsData <- function(bugdir, name, fracs, UmisPerCellLimit = 200, fig_data_path = figure_data_path) {
  #Generate data
  print(paste0("Generating data for ", name))
  bug = readBug(bugdir)
  #filter out low quality cells
  #Should have more than 200 UMIs
  UMIsPerCell = bug %>% group_by(barcode) %>% tally()
  sum(UMIsPerCell$n > UmisPerCellLimit)
  filtBug = bug[bug$barcode %in% UMIsPerCell$barcode[UMIsPerCell$n > UmisPerCellLimit],]

  #skip mitochondrial content for now; don't think it matters for this application
  #should have mitochondrial content of less than 10%
  #genelist = filtBug %>% group_by(barcode) %>% do(genes=c(.$gene))

  #mc = rep(0,length(genelist$barcode))

  #for (i in 1:length(genelist)) {
  #  tmp <- grep(pattern = "^mt-", x = genelist$genes[[i]], value = TRUE)
  #  mc[i] = length(tmp) /length(genelist$genes[[i]])
  #}


  dsBugs = downsampleBugs(filtBug, fracs)
  stats = getStatsFromBugs(dsBugs, fracs, name)

  #create folder
  subpath = paste0(fig_data_path, name, "/")
  if (!file.exists(paste0(fig_data_path, name))){ #there must be no slash at the end here...
    dir.create(file.path(fig_data_path, name))
  }

  #generate proper variable names before saving
  fracsName = paste0("fracs",name)
  bugName = paste0("bug",name)
  bugsName = paste0("dsBugs",name)
  statsName = paste0("stats",name)

  print("Saving BUG...")
  namedSave(list(filtBug, fracs), list(bugName, fracsName), paste0(subpath, "Bug.RData"))
  print("Saving BUGs...")
  namedSave(list(dsBugs, fracs), list(bugsName, fracsName), paste0(subpath, "DsBugs.RData"))
  #save(get(bugsName), get(fracsName), file=paste0(subpath, "DsBugs.RData"))
  print("Saving stats...")
  namedSave(list(stats), list(statsName), paste0(subpath, "Stats.RData"))
  #save(list, file=paste0(subpath, "Counts.RData"))
  #save(get(cpmsName), file=paste0(subpath, "Cpms.RData"))
}

