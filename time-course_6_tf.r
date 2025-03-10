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

#load in the data
files <- c("day1r1.tabular", "day1r2.tabular", "day1r3.tabular", 
           "day3r1.tabular", "day3r2.tabular", "day3r3.tabular", 
           "day6r1.tabular", "day6r2.tabular", "day6r3.tabular", 
           "day9r2.tabular", "day9r3.tabular")
#have a look
read.delim(files[1], nrow=5)

#put into a matrix
x <- readDGE(files, columns=c(3,2))
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

#filter out low-expressed genes
#filtering by 10 reads, going from 23169 to 15858 genes
nrow(deseqset)
keep <- rowSums(counts(deseqset) >= 10) >= 3
deseqsetf <-deseqset[keep,]
nrow(deseqsetf)



## genelist
# using 6 qpcr genes
q_6 <- readLines("6_qpcr.txt")
head(q_6)


# Loop through genes (or create a custom dataframe)
ages_list <- lapply(q_6, function(gene) {
  plotCounts(deseqsetf, gene, intgroup = "condition", returnData = TRUE)
})

# Combine all gene data frames into one
ages_combined <- do.call(rbind, ages_list)
ages_combined$gene <- rep(q_6, each = nrow(ages_list[[1]]))  # Add gene names as a new column


# Plot with ggplot and facet by gene
ggplot(ages_combined, aes(x = as.numeric(as.character(condition)), y = count, colour = gene, group = gene)) + 
  stat_summary(fun = "mean", geom = "line", size = 1) +
  geom_point() +
  labs(y = "Normalized counts (log)", x = "Age (days)") +
  theme_classic() +
  facet_wrap(~gene) +  # This will create separate plots for each gene
  scale_x_continuous(breaks = c(1,3,6,9)) +
  scale_y_log10()



