# Set working directory

# import libraries we will need in this session
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)

getwd()

# Import counts matrix
countsTable <- read.table("DE_counts_F3.txt", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)
# [1] 25279    16

countsTableRound <- round(countsTable) #because DESeq2 doesn't like decimals (and Salmon outputs decimals)
head(countsTableRound)
# table shows how many reads mapped to a certain transcript
# samples on top, transcripts on left side

# import the sample description table
conds <- read.delim("RT_tonsa_F3_samples.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
head(conds)

# How many reads do we have from each sample?
colSums(countsTableRound)
mean(colSums(countsTableRound))
# [1] 17727744
barplot(colSums(countsTableRound), names.arg=colnames(countsTableRound), cex.names=0.5, las=3, ylim=c(0,20000000))
abline(h=mean(colSums(countsTableRound)), col="blue",lwd=2)

# The average number of counts per gene
rowSums(countsTableRound)
mean(rowSums(countsTableRound)) # [1] 11220.54
median(rowSums(countsTableRound)) # [1] 2144
# We see that the mean is much higher than the median. This is driven by a few very highly expressed genes

apply(countsTableRound,2,mean) # 2 in the apply function does the action across columns
apply(countsTableRound,1,mean) # 1 in the apply function does the action across rows

# make histogram for average number of reads per gene
hist(apply(countsTableRound,1,mean), xlim=c(0,1000), breaks=10000)
# many many genes have a very low amount of reads

# create a DESeq object and define the experimental design here with the tilda
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, design= ~ line + environment + line:environment) 

dim(dds)
# same number of genes and samples as counts table

# Filter out genes with too few reads - keep reads with average > 10 reads per sample
dds <- dds[rowSums(counts(dds)) > 160]
dim(dds)
# still same number

# Run the DESeq model to test for differential gene expression
dds <- DESeq(dds)

# List the results you've generated
resultsNames(dds)
#[1] "Intercept"                 
# [2] "line_combined_vs_ambient"  
# [3] "environment_HH_vs_AA"      
# [4] "linecombined.environmentHH"

# PCA to visualize global gene expression patterns
vsd <- vst(dds, blind=FALSE)
data <- plotPCA(vsd, intgroup=c("line","environment"),returnData=TRUE)
percentVar <- round(100 * attr(data,"percentVar"))

# Plot this PCA using ggplot
ggplot(data, aes(PC1,PC2, color=environment, shape=line)) + geom_point(size=4, alpha=0.85) + xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()
# It seems like environment has more of an effect on genetic variance than line. However, ambient is grouped together more than combined line. 

### Order and summarize results from specific contrasts
resInteraction <- results(dds, alpha=0.05)
resInteraction <- resInteraction[order(resInteraction$padj),]
head(resInteraction)
# log2 fold change (MLE): linecombined.environmentHH 
# Wald test p-value: linecombined.environmentHH 
# DataFrame with 6 rows and 6 columns
# baseMean
# <numeric>
#   TRINITY_DN142181_c0_g4    384.448
# TRINITY_DN143012_c0_g4    467.286
# TRINITY_DN131723_c0_g1   1926.655
# TRINITY_DN142181_c0_g18   364.882
# TRINITY_DN145818_c5_g1    297.741
# TRINITY_DN135177_c0_g1   5854.210
# log2FoldChange
# <numeric>
#   TRINITY_DN142181_c0_g4         3.11803
# TRINITY_DN143012_c0_g4        -3.31564
# TRINITY_DN131723_c0_g1         2.60636
# TRINITY_DN142181_c0_g18        3.11752
# TRINITY_DN145818_c5_g1         1.89854
# TRINITY_DN135177_c0_g1         2.38301
# lfcSE
# <numeric>
#   TRINITY_DN142181_c0_g4   0.398309
# TRINITY_DN143012_c0_g4   0.475434
# TRINITY_DN131723_c0_g1   0.387730
# TRINITY_DN142181_c0_g18  0.468239
# TRINITY_DN145818_c5_g1   0.295160
# TRINITY_DN135177_c0_g1   0.396707
# stat
# <numeric>
#   TRINITY_DN142181_c0_g4    7.82817
# TRINITY_DN143012_c0_g4   -6.97391
# TRINITY_DN131723_c0_g1    6.72209
# TRINITY_DN142181_c0_g18   6.65797
# TRINITY_DN145818_c5_g1    6.43225
# TRINITY_DN135177_c0_g1    6.00699
# pvalue
# <numeric>
#   TRINITY_DN142181_c0_g4  4.95009e-15
# TRINITY_DN143012_c0_g4  3.08244e-12
# TRINITY_DN131723_c0_g1  1.79140e-11
# TRINITY_DN142181_c0_g18 2.77627e-11
# TRINITY_DN145818_c5_g1  1.25730e-10
# TRINITY_DN135177_c0_g1  1.89001e-09
# padj
# <numeric>
#   TRINITY_DN142181_c0_g4  1.12976e-10
# TRINITY_DN143012_c0_g4  3.51753e-08
# TRINITY_DN131723_c0_g1  1.36283e-07
# TRINITY_DN142181_c0_g18 1.58407e-07
# TRINITY_DN145818_c5_g1  5.73907e-07
# TRINITY_DN135177_c0_g1  7.18927e-06

### summary of this interaction
summary(resInteraction)
# out of 25279 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 271, 1.1%
# LFC < 0 (down)     : 60, 0.24%
# outliers [1]       : 5, 0.02%
# low counts [2]     : 2451, 9.7%

######### TEST FOR EFFECT OF ENVIRONMENT
#######################
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design = ~ line + environment)

dds <- DESeq(dds, test="LRT", reduced=~line)

# List the results you've generated
resultsNames(dds)
# [1] "Intercept"               
# [2] "line_combined_vs_ambient"
# [3] "environment_HH_vs_AA" 

# Order and list and summarize results from specific contrasts
resEnv <- results(dds, alpha = 0.05)
resEnv <- resEnv[order(resEnv$padj),]
head(resEnv)
# log2 fold change (MLE): environment HH vs AA 
# LRT p-value: '~ line + environment' vs '~ line' 
# DataFrame with 6 rows and 6 columns
# baseMean
# <numeric>
#   TRINITY_DN121599_c1_g1     728.155
# TRINITY_DN150588_c1_g2    1337.692
# TRINITY_DN136932_c13_g13   132.738
# TRINITY_DN146851_c0_g5     465.046
# TRINITY_DN145745_c0_g5     130.345
# TRINITY_DN83766_c0_g1      109.739
# log2FoldChange
# <numeric>
#   TRINITY_DN121599_c1_g1         -6.02233
# TRINITY_DN150588_c1_g2          2.09424
# TRINITY_DN136932_c13_g13        2.09597
# TRINITY_DN146851_c0_g5          1.33330
# TRINITY_DN145745_c0_g5          1.45951
# TRINITY_DN83766_c0_g1          -1.57117
# lfcSE
# <numeric>
#   TRINITY_DN121599_c1_g1    0.472156
# TRINITY_DN150588_c1_g2    0.271553
# TRINITY_DN136932_c13_g13  0.291509
# TRINITY_DN146851_c0_g5    0.193383
# TRINITY_DN145745_c0_g5    0.217831
# TRINITY_DN83766_c0_g1     0.233503
# stat
# <numeric>
#   TRINITY_DN121599_c1_g1    109.3032
# TRINITY_DN150588_c1_g2     53.2838
# TRINITY_DN136932_c13_g13   48.0239
# TRINITY_DN146851_c0_g5     46.1114
# TRINITY_DN145745_c0_g5     43.6364
# TRINITY_DN83766_c0_g1      43.4019
# pvalue
# <numeric>
#   TRINITY_DN121599_c1_g1   1.39269e-25
# TRINITY_DN150588_c1_g2   2.88676e-13
# TRINITY_DN136932_c13_g13 4.21047e-12
# TRINITY_DN146851_c0_g5   1.11718e-11
# TRINITY_DN145745_c0_g5   3.95424e-11
# TRINITY_DN83766_c0_g1    4.45756e-11
# padj
# <numeric>
#   TRINITY_DN121599_c1_g1   3.44996e-21
# TRINITY_DN150588_c1_g2   3.57554e-09
# TRINITY_DN136932_c13_g13 3.47673e-08
# TRINITY_DN146851_c0_g5   6.91867e-08
# TRINITY_DN145745_c0_g5   1.84038e-07
# TRINITY_DN83766_c0_g1    1.84038e-07

summary(resEnv)
# out of 25279 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 513, 2%
# LFC < 0 (down)     : 315, 1.2%
# outliers [1]       : 16, 0.063%
# low counts [2]     : 491, 1.9%

## Up = more expressed in HH environment, down = more expressed in AA.

# exlude na values in column
resEnv <- resEnv[!is.na(resEnv$padj),]
dim(resEnv)
# [1] 24772     6

# differentially expressed genes from the environment
degsEnv <- row.names(resEnv[resEnv$padj < 0.05,]) 
length(degsEnv)
# [1] 828
# 828 is sum of genes upregulated and downregulated.



#########  TEST FOR EFFECT OF LINE
#######################

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design = ~ environment + line)

dds <- DESeq(dds, test="LRT", reduced=~environment)
resultsNames(dds)
# [1] "Intercept"               
# [2] "environment_HH_vs_AA"    
# [3] "line_combined_vs_ambient"
resLine <- results(dds, alpha = 0.05)
resLine <- resLine[order(resLine$padj),]
head(resLine)
# log2 fold change (MLE): line combined vs ambient 
# LRT p-value: '~ environment + line' vs '~ environment' 
# DataFrame with 6 rows and 6 columns
# baseMean
# <numeric>
#   TRINITY_DN132194_c0_g1  1229.949
# TRINITY_DN121089_c0_g2   632.855
# TRINITY_DN134798_c0_g1   152.989
# TRINITY_DN129890_c0_g4   153.602
# TRINITY_DN147342_c0_g4   151.496
# TRINITY_DN134960_c1_g9  2869.532
# log2FoldChange
# <numeric>
#   TRINITY_DN132194_c0_g1        1.53901
# TRINITY_DN121089_c0_g2        1.50181
# TRINITY_DN134798_c0_g1       -2.05171
# TRINITY_DN129890_c0_g4       -3.22080
# TRINITY_DN147342_c0_g4        1.86674
# TRINITY_DN134960_c1_g9        1.94038
# lfcSE
# <numeric>
#   TRINITY_DN132194_c0_g1  0.154825
# TRINITY_DN121089_c0_g2  0.163775
# TRINITY_DN134798_c0_g1  0.224453
# TRINITY_DN129890_c0_g4  0.341257
# TRINITY_DN147342_c0_g4  0.237816
# TRINITY_DN134960_c1_g9  0.245758
# stat
# <numeric>
#   TRINITY_DN132194_c0_g1   94.6814
# TRINITY_DN121089_c0_g2   80.8863
# TRINITY_DN134798_c0_g1   77.6336
# TRINITY_DN129890_c0_g4   76.5916
# TRINITY_DN147342_c0_g4   58.6985
# TRINITY_DN134960_c1_g9   57.2414
# pvalue
# <numeric>
#   TRINITY_DN132194_c0_g1 2.23631e-22
# TRINITY_DN121089_c0_g2 2.39089e-19
# TRINITY_DN134798_c0_g1 1.24040e-18
# TRINITY_DN129890_c0_g4 2.10231e-18
# TRINITY_DN147342_c0_g4 1.83783e-14
# TRINITY_DN134960_c1_g9 3.85471e-14
# padj
# <numeric>
#   TRINITY_DN132194_c0_g1 5.43021e-18
# TRINITY_DN121089_c0_g2 2.90279e-15
# TRINITY_DN134798_c0_g1 1.00398e-14
# TRINITY_DN129890_c0_g4 1.27620e-14
# TRINITY_DN147342_c0_g4 8.92522e-11
# TRINITY_DN134960_c1_g9 1.56000e-10

summary(resLine)
# out of 25279 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 808, 3.2%
# LFC < 0 (down)     : 837, 3.3%
# outliers [1]       : 16, 0.063%
# low counts [2]     : 981, 3.9%
# (mean count < 21)

resLine <- resLine[!is.na(resLine$padj),]
degsline <- row.names(resLine[resLine$padj < 0.05,])
length(degsline)
# [1] 1645

########  TEST FOR INTERACTION
#######################
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design = ~ environment + line + environment:line)

dds <- DESeq(dds, test="LRT", reduced=~environment + line)
resultsNames(dds)
# [1] "Intercept"                 
# [2] "environment_HH_vs_AA"      
# [3] "line_combined_vs_ambient"  
# [4] "environmentHH.linecombined"
resInt <- results(dds, alpha = 0.05)
resInt <- resInt[order(resInt$padj),]
head(resInt)
# log2 fold change (MLE): environmentHH.linecombined 
# LRT p-value: '~ environment + line + environment:line' vs '~ environment + line' 
# DataFrame with 6 rows and 6 columns
# baseMean
# <numeric>
#   TRINITY_DN142181_c0_g4    384.448
# TRINITY_DN143012_c0_g4    467.286
# TRINITY_DN131723_c0_g1   1926.655
# TRINITY_DN142181_c0_g18   364.882
# TRINITY_DN145818_c5_g1    297.741
# TRINITY_DN135177_c0_g1   5854.210
# log2FoldChange
# <numeric>
#   TRINITY_DN142181_c0_g4         3.11803
# TRINITY_DN143012_c0_g4        -3.31564
# TRINITY_DN131723_c0_g1         2.60636
# TRINITY_DN142181_c0_g18        3.11752
# TRINITY_DN145818_c5_g1         1.89854
# TRINITY_DN135177_c0_g1         2.38301
# lfcSE
# <numeric>
#   TRINITY_DN142181_c0_g4   0.398309
# TRINITY_DN143012_c0_g4   0.475434
# TRINITY_DN131723_c0_g1   0.387730
# TRINITY_DN142181_c0_g18  0.468239
# TRINITY_DN145818_c5_g1   0.295160
# TRINITY_DN135177_c0_g1   0.396707
# stat
# <numeric>
#   TRINITY_DN142181_c0_g4    59.3373
# TRINITY_DN143012_c0_g4    46.3684
# TRINITY_DN131723_c0_g1    43.7887
# TRINITY_DN142181_c0_g18   42.7430
# TRINITY_DN145818_c5_g1    40.8415
# TRINITY_DN135177_c0_g1    35.1139
# pvalue
# <numeric>
#   TRINITY_DN142181_c0_g4  1.32838e-14
# TRINITY_DN143012_c0_g4  9.79853e-12
# TRINITY_DN131723_c0_g1  3.65817e-11
# TRINITY_DN142181_c0_g18 6.24250e-11
# TRINITY_DN145818_c5_g1  1.65090e-10
# TRINITY_DN135177_c0_g1  3.10977e-09
# padj
# <numeric>
#   TRINITY_DN142181_c0_g4  2.96667e-10
# TRINITY_DN143012_c0_g4  1.09415e-07
# TRINITY_DN131723_c0_g1  2.72326e-07
# TRINITY_DN142181_c0_g18 3.48535e-07
# TRINITY_DN145818_c5_g1  7.37393e-07
# TRINITY_DN135177_c0_g1  1.15751e-05

summary(resInt)
# out of 25279 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 235, 0.93%
# LFC < 0 (down)     : 48, 0.19%
# outliers [1]       : 5, 0.02%
# low counts [2]     : 2941, 12%
# (mean count < 28)

## Less than 2% of genes showed an interaction
resInt <- resInt[!is.na(resInt$padj),]
degsInt <- row.names(resInt[resInt$padj < 0.05,])
length(degsInt)
# [1] 283

### Plot Individual genes ### 
# Counts of specific top interaction gene! (important validation that the normalization, model is working)
d <-plotCounts(dds, gene="TRINITY_DN142181_c0_g4", intgroup = (c("line","environment")), returnData=TRUE)
d

p <-ggplot(d, aes(x=environment, y=count, color=line, shape=line, group=line)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "line")
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
p


### Plot specific top line effect gene!! ###

d <-plotCounts(dds, gene="TRINITY_DN132194_c0_g1", intgroup = (c("line","environment")), returnData=TRUE)
d

p <-ggplot(d, aes(x=environment, y=count, color=line, shape=line, group=line)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "line")
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
p
### Higher expression in combined, lower expression in ambient. 
## There is a main effect of line!

## Plot specific top environment effect gene!

d <-plotCounts(dds, gene="TRINITY_DN121599_c1_g1", intgroup = (c("line","environment")), returnData=TRUE)
d

p <-ggplot(d, aes(x=environment, y=count, color=line, shape=line, group=line)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "line")
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
p
## In ambient conditions, higher expression in ambient, lower in combined. In combined conditions, they are the same. Both counts at zero.

####### Overlapping Differentially expressed genes in Venn Diagram ###########
library(eulerr)

# Total
length(degsEnv)  # 828
length(degsline)  # 1645
length(degsInt)  # 283

# Intersections
length(intersect(degsEnv,degsline)) # 141
length(intersect(degsEnv,degsInt))  # 14
length(intersect(degsInt,degsline))  # 32

intEL <- intersect(degsEnv,degsline)
length(intersect(degsInt,intEL)) # 7 shared between all 3

# Number unique
828-141-14-7 # 666
1645-141-32-7 # 1465
283-14-32-7 # 230

fit1 <- euler(c("Env" = 666, "Line" = 1465, "Interaction" = 230, "Env&Line" = 141, "Env&Interaction" = 14, "Line&Interaction" = 32, "Env&Line&Interaction" = 7))

plot(fit1,  lty = 1:3, quantities = TRUE)

plot(fit1, quantities = TRUE, fill = "transparent",
     lty = 1:3,
     labels = list(font = 4))

# Heatmap of top 20 genes sorted by pvalue

library(pheatmap)

# Interaction

topgenes <- head(rownames(resInt),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("line","environment")])
pheatmap(mat, annotation_col=df)

# By line
## This heatmap shows that if you sort by line you get very distinct regions that are up and downregulated 
topgenes <- head(rownames(resLine),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("line","environment")])
pheatmap(mat, annotation_col=df)

# By environment
## You get pretty distinct grouping for environment as well, but not as extreme colors, showing not as extreme differential expression
topgenes <- head(rownames(resEnv),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("line","environment")])
pheatmap(mat, annotation_col=df)

