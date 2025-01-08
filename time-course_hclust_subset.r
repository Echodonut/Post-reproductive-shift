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

##############
#transform rnaseq reads to correct for a dependency between the mean and the variance

vsd <- vst(deseqset, blind=FALSE)
head(assay(vsd),3)
vsd_vec <- assay(vsd)
trans_cts <- data.frame(gene = row.names(vsd_vec), vsd_vec)


#we now have the corrected reads as the vector "vsd_vec", which is the same as the counts_transformed.csv
sample_info <- read_csv("./hclust_tut_data/sample_info_cbr.csv")


#select only 1-6 significant genes for example
up1_6 <- read_csv("./hclust_tut_data/up1_6_cbr.csv")

up1_6 <- up1_6 %>% 
  pull(gene)            # extract the gene column as a vector

# Summarise counts 
trans_cts_mean <- trans_cts %>% 
  # convert to long format
  pivot_longer(cols = day1r1:day9r3, names_to = "sample", values_to = "cts")  %>% 
  # join with sample info table
  full_join(sample_info, by = ("sample")) %>% 
  # filter to retain only genes of interest
  filter(gene %in% up1_6) %>% 
  # for each gene
  group_by(gene) %>% 
  # scale the cts column
  mutate(cts_scaled = (cts - mean(cts))/sd(cts)) %>% 
  # for each gene and day
  group_by(gene, day) %>%
  # calculate the mean (scaled) cts
  summarise(mean_cts_scaled = mean(cts_scaled),
            nrep = n()) %>% 
  ungroup()
#look at result
head(trans_cts_mean)
dim(trans_cts_mean)

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

#select only candidate genes
hclust_matrix <- hclust_matrix[up1_6, ]

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
plot(gene_hclust, labels = FALSE, 
     main="", 
     sub="", 
     xlab ="Distance between genes")

abline(h = 4.75, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram


#to sort the genes by clusters, pick any height in the tree and cut it there.
#genes will be assorted to their dendrogram branch.
gene_cluster <- cutree(gene_hclust, k=9) %>%
  #make into tibble
  enframe() %>%
  #rename columns
  rename(gene = name, cluster = value)

#visualize
trans_cts_cluster <- trans_cts_mean %>%
  inner_join(gene_cluster, by = "gene")
head(trans_cts_cluster)



#make datafrom of gene numbers, to add to plot
cluster_sizes <- c("Cluster 1" = 2213, "Cluster 2" = 1540, "Cluster 3" = 104, 
                                      "Cluster 4" = 29, "Cluster 5" = 621, "Cluster 6" = 241,
                                      "Cluster 7" = 67, "Cluster 8" = 201, "Cluster 9" = 58)

cluster_sizes_df <- data.frame(
  cluster = names(cluster_sizes),
  gene_count = as.numeric(cluster_sizes))

# Merge the cluster sizes data into your original dataframe
trans_cts_cluster$cluster <- as.character(trans_cts_cluster$cluster)
trans_cts_cluster <- trans_cts_cluster %>%
  left_join(cluster_sizes_df, by = "cluster")

# Create a labeller function that combines the cluster name with gene count
cluster_labeller <- function(variable, value) {
  # Create a named vector that includes the cluster name and gene count
  labels <- paste("Cluster", value, ": ", cluster_sizes_df$gene_count[match(value, cluster_sizes_df$cluster)], " genes")
  return(labels)
}

#write genes into file. un-comment if you need this.
#clust_count <- gene_cluster %>% count(cluster)
#write.csv(clust_count, "./hclust_tut_data/clust_count_11.csv")

#the plot
#I made alpha value very low, to show by transparency where most lines overlap
p <- trans_cts_cluster %>% 
  ggplot(aes(day, mean_cts_scaled)) +
  geom_line(aes(group = gene), 
            alpha = 0.1,
            linewidth = 0.2,
  ) +
  geom_line(stat = "summary", 
            fun = "median", 
            colour = "brown", 
            size = 1.5, 
            aes(group = 1)
  ) +
  scale_x_continuous(breaks = c(1,3,6,9)) +
  labs(y = "Scaled counts", x = "Age (days)") +
  theme_classic() +
  facet_grid(cols = vars(cluster)) +
theme(
  strip.text = element_text(
    size = 12, 
    face = "bold", 
    hjust = 0.5, 
    vjust = 2.5  # Adjust vertical position of the text above the facet
  )
) 

print(p)