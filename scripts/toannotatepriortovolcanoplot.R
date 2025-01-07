# set working directory
#setwd()

#read 1dayCS file
day1_CS <- read.csv("1dayCSvsControl_differential_expression_allresults.csv", row.names = 1)
head(day1_CS)

# Add a new column 'Diffexpressed' based on conditions
day1_CS$Diffexpressed <- ifelse(day1_CS$padj <= 0.05 & day1_CS$log2FoldChange > 0.585, "Upregulated",
                           ifelse(day1_CS$padj <= 0.05 & day1_CS$log2FoldChange < -0.585, "Downregulated", "Not Significant"))

# Add a new column 'neglogpadj' to calculate -log10(padj)
day1_CS$neglogpadj <- -log10(day1_CS$padj)

View(day1_CS)

# Save results
write.csv(day1_CS, file = "1dayCSvsControl_differential_expression_allresults_forvolcano.csv")

#This file has NAs is pvalue or padj columns. Read why here: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA

# Omit rows with NA in padj column
day1CSnew<- day1_CS[!is.na(day1_CS$padj), ]

# Save new results
write.csv(day1CSnew, file = "1dayCSvsControl_differential_expression_allresults_forvolcanowithoutNA.csv")

#as dataframe
day1CSnew <- as.data.frame(day1CSnew)

#read gene description and ensembl file
names_LTandST <- read.csv("GSE76205_ensemblidonly_MEedited.csv", row.names = 1)
head(names_LTandST)

#Ensure that the Ensembl IDs in day1CSnew are correctly set as row names:
day1CSnew$ensembl_id <- rownames(day1CSnew)

#Ensure that the Ensembl IDs in names are correctly set as row names:
names_LTandST $ensembl_id <- rownames(names_LTandST)

# Merge based on 'ensembl_id'
merged_df <- merge(day1CSnew, names_LTandST , by = "ensembl_id", all.x = TRUE)

# Save results
write.csv(merged_df, file = "Final_1dayCSvsControl_differential_expression_allresults_forvolcano.csv")


#####################################################################################################################################

#read 7dayCS file
day7_CS <- read.csv("7dayCSvsControl_differential_expression_allresults.csv", row.names = 1)
head(day7_CS)

# Add a new column 'Diffexpressed' based on conditions
day7_CS$Diffexpressed <- ifelse(day7_CS$padj <= 0.05 & day7_CS$log2FoldChange > 0.585, "Upregulated",
                                ifelse(day7_CS$padj <= 0.05 & day7_CS$log2FoldChange < -0.585, "Downregulated", "Not Significant"))

# Add a new column 'neglogpadj' to calculate -log10(padj)
day7_CS$neglogpadj <- -log10(day7_CS$padj)

View(day7_CS)

# Save results
write.csv(day7_CS, file = "7dayCSvsControl_differential_expression_allresults_forvolcano.csv")

#This file has NAs is pvalue or padj columns. Read why here: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA

# Omit rows with NA in padj column
day7CSnew<- day7_CS[!is.na(day7_CS$padj), ]

# Save new results
write.csv(day7CSnew, file = "7dayCSvsControl_differential_expression_allresults_forvolcanowithoutNA.csv")

#as dataframe
day7CSnew <- as.data.frame(day7CSnew)

#read gene description and ensembl file
names_LTandST <- read.csv("GSE76205_ensemblidonly_MEedited.csv", row.names = 1)
head(names_LTandST)

#Ensure that the Ensembl IDs in day7CSnew are correctly set as row names:
day7CSnew$ensembl_id <- rownames(day7CSnew)

#Ensure that the Ensembl IDs in names are correctly set as row names:
names_LTandST $ensembl_id <- rownames(names_LTandST)

# Merge based on 'ensembl_id'
merged_df <- merge(day7CSnew, names_LTandST , by = "ensembl_id", all.x = TRUE)

# Save results
write.csv(merged_df, file = "Final_7dayCSvsControl_differential_expression_allresults_forvolcano.csv")


#####################################################################################################################################

#read 1monthCS file
mon1_CS <- read.csv("1monthCSvsControl_differential_expression_allresults.csv", row.names = 1)
head(mon1_CS)

# Add a new column 'Diffexpressed' based on conditions
mon1_CS$Diffexpressed <- ifelse(mon1_CS$padj <= 0.05 & mon1_CS$log2FoldChange > 0.585, "Upregulated",
                                ifelse(mon1_CS$padj <= 0.05 & mon1_CS$log2FoldChange < -0.585, "Downregulated", "Not Significant"))

# Add a new column 'neglogpadj' to calculate -log10(padj)
mon1_CS$neglogpadj <- -log10(mon1_CS$padj)

View(mon1_CS)

# Save results
write.csv(mon1_CS, file = "1monthCSvsControl_differential_expression_allresults_forvolcano.csv")

#This file has NAs is pvalue or padj columns. Read why here: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA

# Omit rows with NA in padj column
mon1CSnew<- mon1_CS[!is.na(mon1_CS$padj), ]

# Save new results
write.csv(mon1CSnew, file = "1monthCSvsControl_differential_expression_allresults_forvolcanowithoutNA.csv")

#as dataframe
mon1CSnew <- as.data.frame(mon1CSnew)

#read gene description and ensembl file
names_LTandST <- read.csv("GSE76205_ensemblidonly_MEedited.csv", row.names = 1)
head(names_LTandST)

#Ensure that the Ensembl IDs in mon1CSnew are correctly set as row names:
mon1CSnew$ensembl_id <- rownames(mon1CSnew)

#Ensure that the Ensembl IDs in names are correctly set as row names:
names_LTandST $ensembl_id <- rownames(names_LTandST)

# Merge based on 'ensembl_id'
merged_df <- merge(mon1CSnew, names_LTandST , by = "ensembl_id", all.x = TRUE)

# Save results
write.csv(merged_df, file = "Final_1monthCSvsControl_differential_expression_allresults_forvolcano.csv")

#####################################################################################################################################

#read 3monthCS file
mon3_CS <- read.csv("3monthCSvsControl_differential_expression_allresults.csv", row.names = 1)
head(mon3_CS)

# Add a new column 'Diffexpressed' based on conditions
mon3_CS$Diffexpressed <- ifelse(mon3_CS$padj <= 0.05 & mon3_CS$log2FoldChange > 0.585, "Upregulated",
                                ifelse(mon3_CS$padj <= 0.05 & mon3_CS$log2FoldChange < -0.585, "Downregulated", "Not Significant"))

# Add a new column 'neglogpadj' to calculate -log10(padj)
mon3_CS$neglogpadj <- -log10(mon3_CS$padj)

View(mon3_CS)

# Save results
write.csv(mon3_CS, file = "3monthCSvsControl_differential_expression_allresults_forvolcano.csv")

#This file has NAs is pvalue or padj columns. Read why here: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA

# Omit rows with NA in padj column
mon3CSnew<- mon3_CS[!is.na(mon3_CS$padj), ]

# Save new results
write.csv(mon3CSnew, file = "3monthCSvsControl_differential_expression_allresults_forvolcanowithoutNA.csv")

#as dataframe
mon3CSnew <- as.data.frame(mon3CSnew)

#read gene description and ensembl file
names_LTandST <- read.csv("GSE76205_ensemblidonly_MEedited.csv", row.names = 1)
head(names_LTandST)

#Ensure that the Ensembl IDs in mon3CSnew are correctly set as row names:
mon3CSnew$ensembl_id <- rownames(mon3CSnew)

#Ensure that the Ensembl IDs in names are correctly set as row names:
names_LTandST $ensembl_id <- rownames(names_LTandST)

# Merge based on 'ensembl_id'
merged_df <- merge(mon3CSnew, names_LTandST , by = "ensembl_id", all.x = TRUE)

# Save results
write.csv(merged_df, file = "Final_3monthCSvsControl_differential_expression_allresults_forvolcano.csv")
#####################################################################################################################################

#read 6monthCS file
mon6_CS <- read.csv("6monthCSvsControl_differential_expression_allresults.csv", row.names = 1)
head(mon6_CS)

# Add a new column 'Diffexpressed' based on conditions
mon6_CS$Diffexpressed <- ifelse(mon6_CS$padj <= 0.05 & mon6_CS$log2FoldChange > 0.585, "Upregulated",
                                ifelse(mon6_CS$padj <= 0.05 & mon6_CS$log2FoldChange < -0.585, "Downregulated", "Not Significant"))

# Add a new column 'neglogpadj' to calculate -log10(padj)
mon6_CS$neglogpadj <- -log10(mon6_CS$padj)

View(mon6_CS)

# Save results
write.csv(mon6_CS, file = "6monthCSvsControl_differential_expression_allresults_forvolcano.csv")

#This file has NAs is pvalue or padj columns. Read why here: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA

# Omit rows with NA in padj column
mon6CSnew<- mon6_CS[!is.na(mon6_CS$padj), ]

# Save new results
write.csv(mon6CSnew, file = "6monthCSvsControl_differential_expression_allresults_forvolcanowithoutNA.csv")

#as dataframe
mon6CSnew <- as.data.frame(mon6CSnew)

#read gene description and ensembl file
names_LTandST <- read.csv("GSE76205_ensemblidonly_MEedited.csv", row.names = 1)
head(names_LTandST)

#Ensure that the Ensembl IDs in mon6CSnew are correctly set as row names:
mon6CSnew$ensembl_id <- rownames(mon6CSnew)

#Ensure that the Ensembl IDs in names are correctly set as row names:
names_LTandST $ensembl_id <- rownames(names_LTandST)

# Merge based on 'ensembl_id'
merged_df <- merge(mon6CSnew, names_LTandST , by = "ensembl_id", all.x = TRUE)

# Save results
write.csv(merged_df, file = "Final_6monthCSvsControl_differential_expression_allresults_forvolcano.csv")
#####################################################################################################################################

#read 9monthCS file
mon9_CS <- read.csv("9monthCSvsControl_differential_expression_allresults.csv", row.names = 1)
head(mon9_CS)

# Add a new column 'Diffexpressed' based on conditions
mon9_CS$Diffexpressed <- ifelse(mon9_CS$padj <= 0.05 & mon9_CS$log2FoldChange > 0.585, "Upregulated",
                                ifelse(mon9_CS$padj <= 0.05 & mon9_CS$log2FoldChange < -0.585, "Downregulated", "Not Significant"))

# Add a new column 'neglogpadj' to calculate -log10(padj)
mon9_CS$neglogpadj <- -log10(mon9_CS$padj)

View(mon9_CS)

# Save results
write.csv(mon9_CS, file = "9monthCSvsControl_differential_expression_allresults_forvolcano.csv")

#This file has NAs is pvalue or padj columns. Read why here: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA

# Omit rows with NA in padj column
mon9CSnew<- mon9_CS[!is.na(mon9_CS$padj), ]

# Save new results
write.csv(mon9CSnew, file = "9monthCSvsControl_differential_expression_allresults_forvolcanowithoutNA.csv")

#as dataframe
mon9CSnew <- as.data.frame(mon9CSnew)

#read gene description and ensembl file
names_LTandST <- read.csv("GSE76205_ensemblidonly_MEedited.csv", row.names = 1)
head(names_LTandST)

#Ensure that the Ensembl IDs in mon9CSnew are correctly set as row names:
mon9CSnew$ensembl_id <- rownames(mon9CSnew)

#Ensure that the Ensembl IDs in names are correctly set as row names:
names_LTandST $ensembl_id <- rownames(names_LTandST)

# Merge based on 'ensembl_id'
merged_df <- merge(mon9CSnew, names_LTandST , by = "ensembl_id", all.x = TRUE)

# Save results
write.csv(merged_df, file = "Final_9monthCSvsControl_differential_expression_allresults_forvolcano.csv")