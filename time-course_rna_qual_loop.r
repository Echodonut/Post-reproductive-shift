setwd("C:/Users/echod/Documents/R_figures/briggsae_Deseq1_featurecounts_data")

library(tximport)
library(tximeta)
library(airway)
library(limma)
library(Glimma)
library(BiocManager)
library(edgeR)
library(DESeq2)
library(ggplot2)
library(tidyverse)

#load in the data
files <- c("day1r1.tabular", "day1r2.tabular", "day1r3.tabular", 
           "day3r1.tabular", "day3r2.tabular", "day3r3.tabular", 
           "day6r1.tabular", "day6r2.tabular", "day6r3.tabular", 
           "day9r2.tabular", "day9r3.tabular")
#have a look
read.delim(files[1], nrow=5)

#put into a matrix
x <- readDGE(files, columns=c(1,2))
class(x)
dim(x)


#rename samples based on their filenames
samplenames <- substring(colnames(x), 1, nchar(colnames(x)))
samplenames

#putting the samples in the matrix x$samples
colnames(x) <- samplenames
group <- as.factor(c("d1", "d1", "d1", "d3", "d3", "d3",
                     "d6", "d6", "d6", "d9", "d9"))
x$samples$group <- group

# make a regular matrix, for use in deseq2. 
countdata <- x [["counts"]]
head(countdata,3)

age = c(1,1,1,3,3,3,6,6,6,9,9)
deseqset <- DESeqDataSetFromMatrix(countData = countdata, 
                                   design = ~ 1, 
                                   colData=data.frame(condition = as.factor(age)))


## genelists. 
#for the loop
looplist <- c("all", "mRNA_processing", "mRNA_splicing", "RNA_localization", "RNA_processing", "RNA_surveillance", "RNA_transport")
i <- 1

#while loop, which selects the various lists of genes and plots their graphs.
#image names are based on the file names
while (i <= length(looplist)) {
      loopfilename <- paste("./RNA_qual/loop/",looplist[i],".png", sep="")
      cat(loopfilename, "\n")
      
      genelistname <- paste("./RNA_qual/RNA_p_genes_csv_",looplist[i],".csv", sep="")
      #or load in different files: mRNA processing, splicing, RNA localization, processing, surveillance, transport
      genelist <- read_csv(genelistname)
      
# extract the gene column as a vector
genelist2 <- genelist %>% 
  pull(genes)
#genelist2 <- substr(genelist2, 4, nchar(genelist2))
RNA <- unique(genelist2)
head(RNA)



# Loop through genes, create a custom dataframe, and skip genes that were filtered out.
ages_list <- lapply(RNA, function(gene) {
  plotCounts(deseqset, gene, intgroup = "condition", returnData = TRUE)
})

# Combine all gene data frames into one
ages_combined <- do.call(rbind, ages_list)
ages_combined$gene <- rep(RNA, each = nrow(ages_list[[1]]))  # Add gene names as a new column



# Plot with ggplot and facet by gene
ggplot(ages_combined, aes(x = as.numeric(as.character(condition)), y = count, group = gene)) + 
  stat_summary(fun = "mean",
               geom = "line",
               size = 1,
               alpha = 0.2) +
  #geom_point() +
  labs(y = "Normalized counts (log)", x = "Age (days)") +
  theme_classic() +

  #facet_wrap(~gene) +  # This will create separate plots for each gene
  scale_x_continuous(breaks = c(1,3,6,9)) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000,100000), limits = c(0.5,100000))+
  geom_line(stat = "summary", 
            fun = "median", 
            colour = "brown", 
            size = 1.5, 
            aes(group = 1)
            )

#print the plot
ggsave(loopfilename, plot = last_plot(), device = "png", 
       width = 800, height = 800, units = "px" )

#increase increment of while loop
i <- i + 1
}


