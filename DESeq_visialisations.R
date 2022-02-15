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
