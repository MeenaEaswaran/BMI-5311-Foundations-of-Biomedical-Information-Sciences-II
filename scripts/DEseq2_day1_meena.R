# set working directory
#setwd()

#install packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

# load libraries
library(DESeq2)

# Step 1: preparing count data ----------------
 
# read in counts data
counts_data <- read.csv("GSE76205_GeneLevel_Raw_data_MEedited_day1.csv", row.names = 1)
head(counts_data)

# Convert count data to integers
counts_data <- round(counts_data) 
head(counts_data)
 
# read in sample info or metadata
colData <- read.csv("day1_metadata_MEedited_Rinputfile.csv", row.names = 1)
View (colData )
 
# making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData)) #should be TRUE
 
#are they in the same order? #Otherwise DEseq2 will cause an error
all(colnames(counts_data) == rownames(colData)) #should be TRUE
 

# Step 2: construct a DESeqDataSet object ----------
 
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                               colData = colData,
                               design = ~ treatment)
 
dds
 
# pre-filtering: removing rows with low gene counts #recommended or not required
# keeping rows that have at least 10 reads total 
keep <- rowSums(counts(dds)) >= 10
keep
dds <- dds[keep,]
 
dds
 
# set the factor level
dds$treatment <- relevel(dds$treatment, ref = "Control")
 
# NOTE: collapse technical replicates; collapsereplicates is provided Deseq only for technical replicates
 
# Step 3: Run DESeq ----------------------
dds <- DESeq(dds)
res <- results(dds)
 
res
# Explore Results ----------------
 
summary(res)
 
#save raw results
write.csv(res, file = "1dayCSvsControl_differential_expression_allresults.csv")
 
# Filter results for significant genes based on padj and fold change considering FC 1.5 and padj 0.05
significant_results <- subset(res, padj <= 0.05 & abs(log2FoldChange) > 0.585)
 
# Print top significant genes
head(significant_results)
 
# Save results
write.csv(significant_results, file = "1dayCSvsControl_differential_expression_significantresultsonly.csv")
