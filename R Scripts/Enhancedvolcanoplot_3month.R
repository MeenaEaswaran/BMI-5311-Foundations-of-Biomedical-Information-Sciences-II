####Enhanced Volcano plot

#Set working directory
#setwd()

#import data
library(readr)
data <- read.csv("Final_3monthCSvsControl_differential_expression_allresults_forvolcano.csv")
head(data)
library(ggplot2)
library(EnhancedVolcano)

#Modify cut-offs for log2FC and P value; specify title; adjust point and label size; remove grid lines; legend positioning
plot1 <- EnhancedVolcano(data,
                lab = data$gene_symbol,
                x = 'log2FoldChange',
                y = 'padj',
                title = '3 months CS-exposed vs Air-exposed Mouse Lungs',
                xlab = expression(bold(log[2] ~ fold ~ change)),
                ylab = expression (bold(''-''~log[10] ~ padj)),
               # xlim = c(-8, 8),
                #ylim = c(0, 15),
                pCutoff = 5.00E-02,
                FCcutoff = 0.585,
                pointSize = 1.0,
                labSize = 3.0,
                labCol = 'black',
               boxedLabels = TRUE,
               colAlpha = 4/5,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                             legendPosition = 'top',
                legendLabSize = 10,
                legendIconSize = 4.0,
               drawConnectors = TRUE,
               widthConnectors = 0.5,
               colConnectors = 'black')


plot1

             
# Save the plot
ggsave("3monthCSvsAir.png", plot = plot1, units = "in", width = 15, height = 10, dpi = 300)
