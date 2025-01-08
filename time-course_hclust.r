
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

##############
#transform rnaseq reads to correct for a dependency between the mean and the variance

vsd <- vst(deseqsetf, blind=FALSE)
head(assay(vsd),3)
vsd_vec <- assay(vsd)

#we now have the corrected reads as the vector "vsd_vec", which is the same as the counts_transformed.csv
trans_cts <- read_csv("./hclust_tut_data/counts_transformed_cbr.csv")
sample_info <- read_csv("./hclust_tut_data/sample_info_cbr.csv")


# Summarise counts 
trans_cts_mean <- trans_cts %>% 
  # convert to long format
  pivot_longer(cols = day1r1:day9r3, names_to = "sample", values_to = "cts")  %>% 
  # join with sample info table
  full_join(sample_info, by = ("sample")) %>% 
  # for each gene
  group_by(gene) %>% 
  # scale the cts column
  mutate(cts_scaled = (cts - mean(cts))/sd(cts)) %>% 
  # for each gene, strain and minute
  group_by(gene, day) %>%
  # calculate the mean (scaled) cts
  summarise(mean_cts_scaled = mean(cts_scaled),
            nrep = n()) %>% 
  ungroup()
#look at result
head(trans_cts_mean)

# plot the relative expressed values
trans_cts_mean %>%
  ggplot(aes(day, mean_cts_scaled)) +
  geom_line(aes(group = gene), alpha = 0.2) + 
  geom_hline(yintercept = 0, colour = "brown", linetype = "dashed") 


## onwards to clustering
# Create a matrix. this changes nothing substantial, just a different variable type
hclust_matrix <- trans_cts %>% 
  select(-gene) %>% 
  as.matrix()
# assign rownames
rownames(hclust_matrix) <- trans_cts$gene

hclust_matrix <- hclust_matrix %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()

#calculate distance between genes (rows)
gene_dist <- dist(hclust_matrix)

#LET'S GET CLUSTERING
gene_hclust <- hclust(gene_dist, method = "complete")

# The default `plot()` function can be used to produce a simple dendrogram
plot(gene_hclust, labels = FALSE)
abline(h = 6.05, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram


#to sort the genes by clusters, pick any height in the tree and cut it there.
#genes will be assorted to their dendrogram branch.
gene_cluster <- cutree(gene_hclust, k=5) %>%
  #make into tibble
  enframe() %>%
  #rename columns
  rename(gene = name, cluster = value)

#visualize
trans_cts_cluster <- trans_cts_mean %>%
  inner_join(gene_cluster, by = "gene")
head(trans_cts_cluster)

#the plot
#I made alpha value very low, to show by transparency where most lines overlap
trans_cts_cluster %>% 
  ggplot(aes(day, mean_cts_scaled)) +
  geom_line(aes(group = gene), 
            alpha = 0.5,
            linewidth = 1,
  ) +
  geom_line(stat = "summary", 
            fun = "median", 
            colour = "brown", 
            size = 1.5, 
            aes(group = 1)
  ) +
  scale_x_continuous(breaks = c(1,3,6,9)) +
  facet_grid(cols = vars(cluster))
