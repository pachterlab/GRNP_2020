#
# Generates the plots in Fig. 1
#
sourcePath = "C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/"

source(paste0(sourcePath,"ButterflyHelpers.R"))

library(ggplot2)
library(ggpubr)

loadStats("EVAL")

#Fig 1B, I - histograms per gene 

colors = gg_color_hue(2)

#read precalculated histograms
h1 = readRDS(paste0(figure_data_path, "Fig1_h1.RDS"))
h2 = readRDS(paste0(figure_data_path, "Fig1_h2.RDS"))

load(file=paste0(figure_data_path, "EVAL/bugEVAL025.RData"))
UMIs25 = bugEVAL025 %>% group_by(gene) %>% tally()
UMIs25$CPM = UMIs25$n*10^6/sum(UMIs25$n)


df1 = data.frame(x = 1:20, y = h1$density[1:20])#30 for full
fig1B_I_1 = ggplot(df1,aes(x=x,y=y)) +
  geom_bar(stat="identity", fill = colors[1]) +
  labs(y="Density", x="Counts per UMI", title="") +
  theme(panel.background = element_rect("white", "white", 0, 
                                          0, "white"))
fig1B_I_1

h2 = hist(collapsedNonFilt$countslist[collapsedNonFilt$gene == "Ubb"][[1]], breaks=seq(0.5, 100.5, by=1), plot=F)
df2 = data.frame(x = 1:20, y = h2$density[1:20])#30 for full
fig1B_I_2 = ggplot(df2,aes(x=x,y=y)) +
  geom_bar(stat="identity", fill = colors[2]) +
  labs(y="Density", x="Counts per UMI", title="") +
  theme(panel.background = element_rect("white", "white", 0, 
                                        0, "white"))
fig1B_I_2

ggsave(
  paste0(figure_path, "fig1B_I_1.png"),
  plot = fig1B_I_1, device = "png",
  width = 3, height = 2.4, dpi = 300)

ggsave(
  paste0(figure_path, "fig1B_I_2.png"),
  plot = fig1B_I_2, device = "png",
  width = 3, height = 2.4, dpi = 300)

#Also generate the CU histograms for Fig 1A

h3 = c(3,1,0,0,0)
df1_1 = data.frame(x = 1:5, y = h3)
fig1A_1 = ggplot(df1_1,aes(x=x,y=y)) +
  geom_bar(stat="identity", fill = colors[1]) +
  scale_y_continuous(breaks = c(0,1,2,3)) +
  labs(y="UMIs", x="Counts per UMI", title="") +
  theme(panel.background = element_rect("white", "white", 0, 
                                        0, "white"))
fig1A_1

ggsave(
  paste0(figure_path, "fig1A_1.png"),
  plot = fig1A_1, device = "png",
  width = 2, height = 2, dpi = 300)


h4 = c(0,2,2,1,0)
df1_2 = data.frame(x = 1:5, y = h4)
fig1A_2 = ggplot(df1_2, aes(x=x,y=y)) +
  geom_bar(stat="identity", fill = colors[2]) +
  scale_y_continuous(breaks = c(0,1,2)) +
  labs(y="UMIs", x="Counts per UMI", title="") +
  theme(panel.background = element_rect("white", "white", 0, 
                                        0, "white"))
fig1A_2

ggsave(
  paste0(figure_path, "fig1A_2.png"),
  plot = fig1A_2, device = "png",
  width = 2, height = 2, dpi = 300)

#predict h3 and h4
dd1_1 = as.matrix(data.frame(1:5,h3));
rSAC1_1 = mod.ztnb.rSAC(dd1_1, incTol = 1e-5, iterIncTol = 200);

rSAC1_1(2)#6.514601


dd1_2 = as.matrix(data.frame(1:5,h4));
rSAC1_2 = mod.ztnb.rSAC(dd1_2, incTol = 1e-5, iterIncTol = 200);
rSAC1_2(2)#5.375271 

#now, fig 1 B II and III

#Fig 1 B II:
############

r1 = data.frame(statsEVAL)[statsEVAL$gene == "Vmn1r13",]
r2 = data.frame(statsEVAL)[statsEVAL$gene == "Ubb",]
r1_25 = UMIs25$CPM[UMIs25$gene == "Vmn1r13"]
r2_25 = UMIs25$CPM[UMIs25$gene == "Ubb"]


data1 = c(r1$CPM_EVAL_d_5, r1$CPM_EVAL_d_10, r1$CPM_EVAL_d_20, r1_25, r1$CPM_EVAL_d_40, r1$CPM_EVAL_d_60, r1$CPM_EVAL_d_80, r1$CPM_EVAL_d_100)
data2 = c(r2$CPM_EVAL_d_5, r2$CPM_EVAL_d_10, r2$CPM_EVAL_d_20, r2_25, r2$CPM_EVAL_d_40, r2$CPM_EVAL_d_60, r2$CPM_EVAL_d_80, r2$CPM_EVAL_d_100)

data=c(data1,data2)
logData = log2(data + 1)
gene = factor(c(rep(0,8),rep(1,8)), c(0,1), c("Vmn1r13", "Ubb"))
x = c(0.05,0.1,0.2,0.25,0.4,0.6,0.8,1);
df = data.frame(data=logData, Gene=gene, x=log2(x))

dfline = data.frame(x=c(log2(0.25), log2(0.25)), data=c(7.8, 9.5) )
dfmark = df[df$x == log2(0.25),]

fig1B_II = ggplot(df, aes(x=x, y=data, colour = Gene)) + 
  geom_line() + geom_point() +
  ylab(expression(Log[2]*"(CPM + 1)")) +
  geom_line(data=dfline, color="black", linetype=2) +
  geom_point(data=dfmark, size=6, shape=1, color="black") +
  #scale_x_continuous(name="Number of reads", breaks = c(log2(0.05), -3, -2, -1, 0), labels = c("1x", "2.5x", "5x", "10x", "20x")) +
  scale_x_continuous(name="Number of reads", breaks = c(log2(0.0625), log2(0.125), log2(0.25), log2(0.5), log2(1)), labels = c("1/16 x", "1/8 x", "1/4 x", "1/2 x", "1 x")) +
  ylim(c(7,10.2)) +
  theme(panel.background = element_rect("white", "white", 0, 
                                        0, "white"), legend.position= "none")

fig1B_II

ggsave(
  paste0(figure_path, "fig1B_II.png"),
  plot = fig1B_II, device = "png",
  width = 3, height = 3, dpi = 300)


#fig1B_III
###############
#create prediction data

#load precalculated prediction data
r1 = readRDS(paste0(figure_data_path, "Fig1_r1_III.RDS"))
r2 = readRDS(paste0(figure_data_path, "Fig1_r2_III.RDS"))

r1_025 = readRDS(paste0(figure_data_path, "Fig1_r1_025_III.RDS"))
r2_025 = readRDS(paste0(figure_data_path, "Fig1_r2_025_III.RDS"))



data1Temp = unlist(c(r1[2:8]))
data2Temp = unlist(c(r2[2:8]))

data1 = c(data1Temp[1:3], r1_025, data1Temp[4:7])
data2 = c(data2Temp[1:3], r2_025, data2Temp[4:7])

data=c(data1,data2)
logData = log2(data + 1)
gene = factor(c(rep(0,8),rep(1,8)), c(0,1), c("Gene 1 - Vmn1r13", "Gene 2 - Ubb"))
x = c(0.05,0.1,0.2,0.25,0.4,0.6,0.8,1);
df = data.frame(data=logData, Gene=gene, x=log2(x))
dfline = data.frame(x=c(log2(0.25), log2(0.25)), data=c(8.3, 10.13) )
dfmark = df[df$x == log2(0.25),]

fig1B_III = ggplot(df, aes(x=x, y=data, colour = Gene)) + 
  geom_line() + geom_point() +
  ylab(expression(Log[2]*"(CPM + 1)")) +
  geom_line(data=dfline, color="black", linetype=2) +
  geom_point(data=dfmark, size=6, shape=1, color="black") +
  ylim(c(7,10.2)) +
  scale_x_continuous(name="Number of reads", breaks = c(log2(0.0625), log2(0.125), log2(0.25), log2(0.5), log2(1)), labels = c("1/16 x", "1/8 x", "1/4 x", "1/2 x", "1 x")) +
  theme(panel.background = element_rect("white", "white", 0, 0, "white"), legend.position= "none")

fig1B_III

ggsave(
  paste0(figure_path, "fig1B_III.png"),
  plot = fig1B_III, device = "png",
  width = 3, height = 3, dpi = 300)



#create a plot with Legend - copy the legend from here...
#uses data from previous plot
fig1B_Legend = ggplot(df, aes(x=x, y=data, colour = Gene)) + 
  geom_line() +
  theme(panel.background = element_rect("white", "white", 0, 
        0, "white"),legend.position= "bottom", legend.direction = "vertical", legend.title = element_blank())

fig1B_Legend

ggsave(
  paste0(figure_path, "fig1B_Legend.png"),
  plot = fig1B_Legend, device = "png",
  width = 3, height = 3, dpi = 300)

#get the data needed for the main text:
dsExpr = statsEVAL$CPM_EVAL_d_5
fullExpr = statsEVAL$CPM_EVAL_d_100
gene1 = statsEVAL$gene == "Vmn1r13"
gene2 = statsEVAL$gene == "Ubb"
fullExpr[gene1]/fullExpr[gene2]#2.404394
dsExpr[gene2]/dsExpr[gene1]#4.058824

#get the full number of reads
length(dsBugsEVAL)
fullBug = dsBugsEVAL[[7]]
dsBug = dsBugsEVAL[[1]]
sum(fullBug$count)#29913038
sum(dsBug$count)#1495651
