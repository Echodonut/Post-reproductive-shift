#R script to generate MA plots of deseq2 for post-repro transcriptome project.
#When running the script for multiple datasets, don't forget to change nums in l16-22 and l43, and fiddle around title position on l90
setwd("C:/Users/echod/Documents/R_figures/briggsae_Deseq1_featurecounts_data")

library(tximport)
library(limma) # essential for drawing MA plots apparently
library(Glimma)
library(BiocManager)
library(edgeR)
library(DESeq2)
library(ggplot2)
library(tidyverse)

#load in the data
files <- c(
 #           "day1r1.tabular", "day1r2.tabular", "day1r3.tabular", 
         "day3r1.tabular", "day3r2.tabular", "day3r3.tabular"
  ,
         "day6r1.tabular", "day6r2.tabular", "day6r3.tabular"
# ,
          )
  #         "day9r2.tabular", "day9r3.tabular")


count_data_list <- lapply(files, function(file) {
  data <- read.table(file, header = TRUE, sep = "\t", row.names = 1)
  data <- data[, 1, drop = FALSE]  # Extract the counts (column 2)
  # Convert the counts to numeric (important step)
  data[] <- lapply(data, as.numeric)
  
  colnames(data) <- sub(".tabular", "", basename(file))  # Remove file extension and name the columns
  return(data)
})

# Combine the data into one data frame (keep gene names as row names)
counts <- do.call(cbind, count_data_list)
# Name the columns as your sample names
# Check the data
head(counts)

sample_info <- data.frame(
  condition = factor(rep(c("Day 3", "Day 6"), times = c(3,3))),
  row.names = colnames(counts)
)
head(sample_info)

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = sample_info, design = ~ condition)

# Check the DESeqDataSet object
dds
dds <- DESeq(dds)

# Get the results (log2 fold changes between the two timepoints)
res <- results(dds)

# View the results
head(res)

setwd("C:/Users/echod/Documents/R_figures/post-repro_revisions_march25/plots")

#write the figure title
fig_tit <- paste((sample_info$condition)[1], "to", (sample_info$condition)[4])
png_tit <- paste ("cbr",fig_tit,".png")
png_tit <- gsub(" ", "", png_tit)

# MA plot to visualize differential expression (Rwindow version, does not look representative)
# if you die in the plot you die in real life
par(mar = c(4.5, 4.5, 1.5, 0.5),font.lab = 2, cex.lab = 1.5, cex.axis =1.5)  # Adjusting the top margin (3rd value) to reduce space between title and plot
plotMA(res, ylim = c(-11, 11))
#putting label in the corner
usr <- par("usr")
text(x = 0.035, y = (usr[4]*0.94), labels = fig_tit, pos = 4, cex = 1.5, font = 2)


# MA plot to visualize differential expression (printed version)
png(png_tit, width = 7, height = 7, units = "in", res = 300)
  par(mar = c(4.5, 4.5, 1.5, 0.5),font.lab = 1, cex.lab = 2, cex.axis =2)  # Adjusting the top margin (3rd value) to reduce space between title and plot
  plotMA(res, ylim = c(-12, 12))

  #putting label in the corner. Need to fiddle around with the values. sorry.
  usr <- par("usr")
  text(x = 0.054, y = (usr[4]*0.95), labels = fig_tit, pos = 4, cex = 2, font = 2)
  
  dev.off()  # Close the device and save the plot

