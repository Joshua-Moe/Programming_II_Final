# BINF-6112-001 Programming II | Spring 2020 | Final Project | Joshua Moe | jmoe5@uncc.edu
# Differential gene expression in mice: sleep deprived v. control   


"First, I am going to set my current working directory to where my project files are located (manually setting it from the menu. Setting it in the command line fails for some reason).
I also will be loading the required packages."


setwd("C:/UNC_Spring_2020/BINF_Programming_II/Projects/SRA_data")

library("pheatmap")
library("DESeq2")
library("ggplot2")
library("calibrate") 





"Now, I will be reading in mm10v2.counts file created in STEPS 1 through 4."

countdata <- read.table("C:/UNC_Spring_2020/BINF_Programming_II/Projects/SRA_data/mm10v2.counts", header = TRUE, row.names = 1)





"I need my countdata to be a matrix, of which I only need the last two columns--columns 6 and 7.
I have to then set my count_matrix as numeric--for whatever reason the 
entire matrix was in string characters."

countdata_matrix <- as.matrix(countdata)
count_matrix = countdata_matrix[,c(6,7)] #This is the matrix that only has the read count columns.
class(count_matrix) <- "numeric"






"Here, I am setting up the prepwork for developing a heatmap that will only display 
the top 25 genes that are above the average overall expression levels."

condition <- factor(c(rep("A",1),rep("B",1)))
colData <- data.frame(row.names = colnames(count_matrix), condition)
dds<- DESeqDataSetFromMatrix(countData = count_matrix, colData = colData, design = ~condition)
select <- order(rowMeans(counts(dds,normalized=FALSE)),decreasing=TRUE)[1:25]






"Now, all that is left is to create the heatmap."

pheatmap(count_matrix[select,])


