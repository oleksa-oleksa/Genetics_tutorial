# Task: DESeq Visualisations
# We start with data loading and preparations
# as it was done in provided tutorial

# Vusialisations can be found after data loading part

readcounts=read.table("Tutorial_macaques_readcounts.txt", header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE)
dim(readcounts)
head(readcounts)

coldata=read.table("Tutorial_macaques_ids.txt", header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE)
dim(coldata)
head(coldata)

# DESeq2 is based on a model using the negative binomial distribution.
library("DESeq2")

# The conditions treatment vs. control (LPS vs. NC) need to be defined as factors for the comparison. 
# The reference (control) in this example is “NC”.
coldata$condition = factor(coldata$condition)

# Make a DESeqDataSet object and look at it. 
# It will hold the count information for each gene for each sample.
dds=DESeqDataSetFromMatrix(countData = readcounts, colData = coldata, design = ~ condition)
dds

# Filter genes with very low expression to reduce noise. 
# Here we will remove all genes with lower than 10 reads in total across all samples.
keep = rowSums(counts(dds)) >= 10
dds = dds[keep,]
dds
# result of filter: rownames(32386) is reduces to rownames(21391)

# To specify the factors for the comparison use the function factor(). In this first example, we will compare NC vs. LPS.
dds$condition = factor(dds$condition, levels = c("NC","LPS"))

# The function DESeq() implements all steps from normalization to the comparison: 
# The first step is always normalization. 
# The core assumption for calculating the dispersion (estimateDispersions) is that the mean is a good predictor 
# of the variance, i.e., that genes with a similar expression level also have similar variance across replicates.
dds = DESeq(dds)

result = results(dds)
result
mcols(result)$description

# To reduce some noise, you can shrink the results. 
resultLFC = lfcShrink(dds, coef="condition_LPS_vs_NC", type="normal")
resultLFC
summary(result)

# count how many genes have an adjusted p-value of smaller than 0.1 or 0.05, respectively.
resultOrdered = result[order(result$pvalue),]
sum(result$padj < 0.1, na.rm=TRUE) # [1] 1607
sum(result$padj < 0.05, na.rm=TRUE) # [1] 1295

# Filtering of genes with False Discovery Rate (FDR) smaller than 0.05 and assigning them to a new object:
result005 = results(dds, alpha=0.05)
summary(result005)

# About MA Plot: In DESeq2, the function plotMA shows the log2 fold changes attributable 
# to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet. 
# Points will be colored red if the adjusted p value is less than 0.1. 
# Points which fall out of the window are plotted as open triangles pointing either up or down.
plotMA(result, ylim=c(-2,2))
# When the shrunken results are plotted with plotMA(), 
# you can see that the noise from lowly expressed genes is reduced.
plotMA(resultLFC, ylim=c(-2,2))

# we can plot the gene with the lowest adjusted p-value. 
# In this plot, counts are normalized by sequencing depth and a pseudocount 
# of 1/2 is added to allow for log scale plotting.
plotCounts(dds, gene=which.min(result$padj), intgroup="condition")

###### 
# Visualisations

# The function plotDispEsts visualizes DESeq2’s dispersion estimates
# The black points are the dispersion estimates for each gene as obtained 
# by considering the informationfrom each gene separately.  
# Unless one has many samples,  these values fluctuate strongly around their true values.
plotDispEsts(dds, ylim = c(1e-6, 1e1))

# Histogram of the pvalues returned by the test for differential expression
hist(resultLFC$pvalue, breaks=20, col="orange" )

hist(result$padj, breaks=40, col="violet" )

###########

# because we are interested in treated vs untreated, we set 'coef=2'
resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resultLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")

#####
# examine the counts of reads for a single gene across the groups. 
# A simple function for making this plot is plotCounts, 
# which normalizes counts by the estimated size factors 
# (or normalization factors if these were used) and adds a pseudocount of 1/2 
# to allow for log scale plotting.
d <- plotCounts(dds, gene=which.min(result$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))


####
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

# These transformation functions return an object of class DESeqTransform which is 
# a subclass of RangedSummarizedExperiment. 
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

meanSdPlot(assay(vsd))


#####
library("pheatmap")
colData(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)["condition"])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
##
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

####
# Sample clustering:  the dist function applied to the transpose 
# of the transformed count matrix to get sample-to-sample distances.
library("RColorBrewer")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# A heatmap of this distance matrix gives us an overview over similarities and dissimilarities 
# between samples. We have to provide a hierarchical clustering hc to the heatmap function 
# based on the sample distances, or else the heatmap function would calculate a clustering 
# based on the distances between the rows/columns of the distance matrix.
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
