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
coldata$condition = relevel(coldata$condition, ref = 'NC') 

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
write.csv(as.data.frame(resultOrdered), file="DESeq2_DEgenes_condition_LPS_NC.csv")

#################################
# edgeR uses empirical Bayes estimation and exact tests based on 
# the negative binomial distribution to call differentially expressed genes.
library("edgeR")

# Make an object with the read count and sample information and look at it. 
# You will see how samples got assigned to groups (LPS or NC), what the library sizes are, the counts for each gene etc.
count_edgeR_obj=DGEList(counts=readcounts, group=coldata$condition)
count_edgeR_obj

# Normalization
count_edgeR_obj=estimateCommonDisp(count_edgeR_obj)
count_edgeR_obj=estimateTagwiseDisp(count_edgeR_obj)

edgeR_DEgenes=exactTest(count_edgeR_obj)
edgeR_DEgenes

# The function topTags() is used to show the top differentially expressed genes (default: based on p-value).
topTags(edgeR_DEgenes)

# To show the top differentially expressed genes based on fold change use:
topTags(edgeR_DEgenes, sort.by = "logFC")

#As seen above, the edgeR_DEgenes object contains multiple elements. 
# The first one is the table with logFC, logCPM, and p-values for each gene. 
# To get access to this table and assign it to a new variable, call:
edgeR_DEgenesTable=edgeR_DEgenes$table
head(edgeR_DEgenesTable)

# extract significant genes
signedgeR_DEgenes=edgeR_DEgenesTable[edgeR_DEgenesTable[,3]<0.05,]

# Write a result file with genes sorted by p-value.
edgeROrdered <- edgeR_DEgenesTable[order(edgeR_DEgenesTable$PValue),]
write.csv(as.data.frame(edgeR_DEgenesTable), file="edgeR_DEgenes_condition_LPS_NC.csv")

# Which method calls more significant genes, DESeq2 or edgeR?
# - edgeR

# colData() shows the different factors of the experimental design. 
colData(dds)

# With unique() you can then ask for a list of all study groups or ranks, respectively.
unique(colData(dds)$study_group)
unique(colData(dds)$rank)
# ranks needs to be converted to factors
coldata$rank = as.factor(coldata$rank)

# We make a new DESeqDataSet object. It will hold the count information for each gene for each sample.
dds_interact=DESeqDataSetFromMatrix(countData = readcounts, 
                                    colData = coldata, 
                                    design = ~ condition + rank + condition:rank)
dds_interact

# filter out genes with very low expression to reduce noise, removing all genes with lower than 10 reads
keep = rowSums(counts(dds_interact)) >= 10
dds_interact = dds_interact[keep,]
dds_interact

# Now we can specify the factors for the comparisons. 
# Here we will consider the condition and the ranks.
dds_interact$condition = factor(dds_interact$condition, levels = c("NC","LPS"))
dds_interact$rank = factor(dds_interact$rank, levels = c("1","2","3","4","5"))

###################################
# Differential expression analysis

dds_interact = DESeq(dds_interact)
results(dds_interact)

# To show the comparisons that were done use the function resultNames() on the DESeq object.
resultsNames(dds_interact)

#########################
require(wTO)
require(magrittr)
NC = readcounts[,coldata$condition == 'NC']
dim(NC)
NC = NC[rowSums(NC)> 10,]
dim(NC)
# Then, collecting all LPS samples and removing genes that have less than 10 counts
LPS = readcounts[,coldata$condition == 'LPS']
dim(LPS)
LPS = LPS[rowSums(LPS)> 10,]
dim(LPS)
# select only the significant genes from the tables with LPS and NC samples. 
# But before that, we will have a look at the results from above again
summary(result)
DE_genes = subset(row.names(result), result$padj<0.01)
NC = subset(NC, row.names(NC) %in% DE_genes)
LPS = subset(LPS, row.names(LPS) %in% DE_genes)
# We want to make a TF wTO network. Thus, we need to retrieve the information about TF genes. 
# The list of TFs is in the file “TFs.txt”. Let’s read this file into our R session.
TFs = read.table("TFs.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
dim(TFs)
TFs = TFs[,1]
length(TFs)

library(biomaRt)
listMarts()
ensembl=useMart("ensembl")

listDatasets(ensembl)
mart = useDataset("mmulatta_gene_ensembl", useMart("ensembl"))

expressedGenes=row.names(result)
listAttributes(mart)

GeneSymbols = getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", 'external_gene_name'),values=expressedGenes,mart= mart)
dim(GeneSymbols)
head (GeneSymbols)

require(plyr) 

LPS$ensembl_gene_id = row.names(LPS)
LPS = join(LPS, GeneSymbols, type = 'inner', match = 'first') 

head(LPS)
# Removing duplicates:
LPS = LPS[!duplicated(LPS$external_gene_name), ]
row.names(LPS) = LPS$external_gene_name
head(LPS)

ncol(LPS)
names(LPS)

# Removing first column with Ensembl IDs and last column with GeneSymbols:
LPS = LPS[, -c(1,27)]
head(LPS)

# Performing the same annotation and reformatting for the NC table:
NC$ensembl_gene_id = row.names(NC)
NC = join(NC, GeneSymbols, type = 'inner', match = 'first')

NC = NC[!duplicated(NC$external_gene_name), ]
row.names(NC) = NC$external_gene_name
head(NC)

names(NC)
