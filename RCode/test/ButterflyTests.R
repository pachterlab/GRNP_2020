#before anything else, you need to setup paths:
#for example
#source("C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/paths.R"))
#or, in colab:
#source("GRNP_2020/RCode/pathsGoogleColab.R")

source(paste0(sourcePath,"ButterflyHelpers.R"))
source(paste0(sourcePath,"preseqHelpers.R"))
source(paste0(sourcePath,"GenBugSummary.R"))
source(paste0(sourcePath,"BUGProcessingHelpers.R"))

test_fig_data_path = paste0(sourcePath, "test/tmp/")

#TCR0001 - ClosestDists:
##############################################
bugTest = read.table(paste0(sourcePath, "test/TestClosestDists.txt"), header = T, stringsAsFactors = F)

all(ClosestDists(bugTest, bugTest[bugTest$gene=="g1",], 10) == c(2,0,1,3,0,0,0,0,0,0)) #ok
all(ClosestDists(bugTest, bugTest[bugTest$gene=="g2",], 10) == c(1,0,0,0,1,0,0,0,0,0)) #ok
all(ClosestDists(bugTest, bugTest[bugTest$gene=="g3",], 10) == c(0,0,0,0,1,2,0,0,0,0)) #ok


#TCR0002 - geneIndices2Symbols
#####################################
#31 == ENSMUSG00000076800.1 == Trav6n-5
res = geneIndices2Symbols(31, paste0(sourcePath, "test/smallBug/bus_output/coll.genes.txt"), paste0(sourcePath, "test/smallBug/bus_output/transcripts_to_genes.txt"))
res == "Trav6n-5" # ok

#TCR0003 - createStandardBugsData
#test 1 simple bug only
#so, there are 8 lines in the bug. The two last lines should be discarded, 
#one because it is multimapped and one because its cell has only one read
#####################################

#create the files for the dataset in the tmp folder
createStandardBugsData(paste0(sourcePath, "test/smallBug/"), "smallBug", c(0.5,1), UmisPerCellLimit = 1, fig_data_path = test_fig_data_path)

#check full bug file briefly
loadBug("smallBug", 1, fig_data_path = test_fig_data_path)
smallBug = getBug("smallBug",1)
all(dim(smallBug) == c(6,4)) #the size of the bug, ok
smallBug[[4,3]] == "Trbd2" #check of one value, ok
rmBug("smallBug", 1)

#check downsampled bug file briefly
loadBug("smallBug", 0.5, fig_data_path = test_fig_data_path)
smallBug50 = getBug("smallBug",0.5)
sum(smallBug50[,4]) == 15 #downsampling, OK
rmBug("smallBug", 0.5)

#check that the stats are ok
loadStats("smallBug", fig_data_path = test_fig_data_path)
res = statssmallBug
res$UMIs_smallBug_d_100[1] == 4   #UMIs Gene 1, ok
res$UMIs_smallBug_d_100[2] == 2   #UMIs Gene 2, ok
res$Counts_smallBug_d_100[1] == 18   #counts Gene 1, ok
res$Counts_smallBug_d_100[2] == 12   #counts Gene 2, ok
abs(res$CPM_smallBug_d_100[1] - 10^6*4/6) < 0.0001   #CPM Gene 1, ok
abs(res$CPM_smallBug_d_100[2] - 10^6*2/6) < 0.0001   #CPM Gene 2, ok
res$FracOnes_smallBug_d_100[1] == 0.25   #fracOnes Gene 1, ok
res$FracOnes_smallBug_d_100[2] == 0   #fracOnes Gene 2, ok
res$CountsPerUMI_smallBug_d_100[1] == 18/4   #counts per UMI Gene 1, ok
res$CountsPerUMI_smallBug_d_100[2] == 12/2   #counts per UMI Gene 2, ok
sum(res$Counts_smallBug_d_50) == 15 #checking downsampling again, ok
res$gene[1] == "Trbd1" #ok
res$gene[2] == "Trbd2" #ok

#TCR0005 - Good-Toulmin
#####################################

dat = c(1,1,1,2,2,3,4)
h = hist(dat, breaks=seq(0.5, 4.5, by=1), plot = F)
goodToulmin(h,2) == 8 #3-2+1-1 + existing number, ok

#TCR0006 - downSampleManyTimesAndGetHist
#####################################

#just check that we get the right total number of counts
loadBug("smallBug", 1, fig_data_path = paste0(sourcePath, "test/tmp/"))
smallBug = getBug("smallBug",1)
histMany = downSampleManyTimesAndGetHist(smallBug, 0.5, numTimes=20)
histAllGenes = colSums(histMany)
totCounts = sum(histAllGenes * 1:100)
totCounts == 30*20*0.5 #The dataset has 30 counts, downsampled to 0.5, repeated 20 times. Ok
rmBug("smallBug", 1)


#TCR0007 - getDsHist
#####################################

#just check that we get the right total number of counts
loadBug("smallBug", 1, fig_data_path = paste0(sourcePath, "test/tmp/"))
smallBug = getBug("smallBug",1)
h = getDsHist(smallBug)
all(h[1,1:10] == c(1,0,0,2,0,0,0,0,1,0)) #ok
all(h[2,1:10] == c(0,0,0,0,0,2,0,0,0,0)) #ok
rmBug("smallBug", 1)

#TCR0008 - downSampleBUGNTimes
#####################################
#test that the expression is halved and multiplied with n if we have only single-copy molecules
bugOnlyOnes = as_tibble(read.table(paste0(sourcePath, "test/bugOnlyOnes.txt"), header=T, stringsAsFactors=F))
m = downSampleBUGNTimes(bugOnlyOnes, 0.5, 5)
sum(m$n) == 25 #ok
#test that if we have one molecule from gene 0 with 3 copies and one from gene 1 with 1 copy, 
#Gene 0 is not affected by downsampling by 50%.
bugUneven = as_tibble(read.table(paste0(sourcePath, "test/bugUneven.txt"), header=T, stringsAsFactors=F))
m = downSampleBUGNTimes(bugUneven, 0.5, 5)
m$n[1] == 5 #ok

#TCR0009 - poolPrediction
#####################################

#We use 3 different genes, which are all the same in the bug (4 UMIs, hist 3 1), not much amplified. 
#Gene1 is well amplified in the pool, and has 100 UMIs
#Gene2 has the same amplification as the bug, and has 100 UMIs
#Gene3 is well amplified in the pool, but only 1 UMI (that is weird, but should work)
#so, Gene1 is not expected to increase much, Gene2 a lot, Gene3 somewhere in between
bugPP = as_tibble(read.table(paste0(sourcePath, "test/bugPooled.txt"), header=T, stringsAsFactors=F))
UMIs = as_tibble(read.table(paste0(sourcePath, "test/pooledUMIs.txt"), header=T, stringsAsFactors=F))
h1 = as.matrix(read.table(paste0(sourcePath, "test/pooledHist1.txt"), header=F, row.names=1, stringsAsFactors=F))
h2 = as.matrix(read.table(paste0(sourcePath, "test/pooledHist2.txt"), header=F, row.names=1, stringsAsFactors=F))
pHList = list(UMIs, list(h1,h2))
pred = poolPrediction(bugPP, t=10, pHList, usePoolLimit = 100000)
pv = pred[[2]]
pv[1] < pv[2]#ok
pv[1] < pv[3]#ok
pv[2] > pv[3]#ok


#TCR0010 - genBugSummary
#####################################
createStandardBugsData(paste0(sourcePath, "test/statsBug/"), "statsBug", c(0.5,1), UmisPerCellLimit = 1, fig_data_path = paste0(sourcePath, "test/tmp/"))
loadBug("statsBug", fig_data_path = paste0(sourcePath, "test/tmp/"))
genBugSummary("statsBug", "Gene1", "Gene2", 10, fig_data_path = paste0(sourcePath, "test/tmp/"))

#now read the summary file and check that it reported the expected values:
con = file(paste0(test_fig_data_path, "statsBug/ds_summary.txt"))
lines = readLines(con)
#print(lines)
close(con)

unlist(strsplit(lines[[3]], "\\s+"))[2] == 13 #num UMIs, ok
unlist(strsplit(lines[[4]], "\\s+"))[2] == 3 #num cells, ok
unlist(strsplit(lines[[5]], "\\s+"))[2] == 37 #counts, ok
unlist(strsplit(lines[[5]], "\\s+"))[2] == 37 #counts, ok
unlist(strsplit(lines[[6]], "\\s+"))[2] == 37/13 #counts per UMI, ok
unlist(strsplit(lines[[7]], "\\s+"))[2] == 13/3 #UMIs per cell, ok
unlist(strsplit(lines[[8]], "\\s+"))[2] == 37/3 #counts per cell, ok
unlist(strsplit(lines[[9]], "\\s+"))[2] == 5/13 #totFracOnes, ok
#Gene1 is low, Gene2 is high

all(unlist(strsplit(lines[[10]], "\\s+"))[2:4] == c("1,","1,","3,")) # f1H
all(unlist(strsplit(lines[[11]], "\\s+"))[2:3] == c("1,","4,")) # f1L
all(unlist(strsplit(lines[[12]], "\\s+"))[2:4] == c("0.2,","0.2,","0.6,")) # f1HFrac
all(unlist(strsplit(lines[[13]], "\\s+"))[2:3] == c("0.2,","0.8,")) # f1LFrac
all(unlist(strsplit(lines[[14]], "\\s+"))[2:4] == c("1,","3,", "1,")) # 1cpy
all(unlist(strsplit(lines[[15]], "\\s+"))[2:4] == c("1,","1,", "1,")) # 2cpy
all(unlist(strsplit(lines[[16]], "\\s+"))[2:4] == c("0,","4,", "1,")) # >3cpy
all(unlist(strsplit(lines[[17]], "\\s+"))[2:4] == c("0.2,","0.6,", "0.2,")) # 1cpy frac
#skip the rest of the frac, it is a trivial calculation and they have a lot of decimals



