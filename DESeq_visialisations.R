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

plotMA(result, ylim=c(-2,2))
# When the shrunken results are plotted with plotMA(), 
# you can see that the noise from lowly expressed genes is reduced.
plotMA(resultLFC, ylim=c(-2,2))

# we can plot the gene with the lowest adjusted p-value. 
# In this plot, counts are normalized by sequencing depth and a pseudocount 
# of 1/2 is added to allow for log scale plotting.
plotCounts(dds, gene=which.min(result$padj), intgroup="condition")
