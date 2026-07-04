#Full R workflow for cutaneous leishmoniasis

# load packages----
library(rhdf5) #provides functions for handling hdf5 file formats (kallisto outputs bootstraps in this format)
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tximport) # package for getting Kallisto results into R
library(ensembldb) #helps deal with ensembl
library(EnsDb.Hsapiens.v86) #v86 is an older version so most of the transcripts cannot be annotated with this package
library(edgeR) # well known package for differential expression analysis, but we only use for the DGEList object and for normalization methods
library(matrixStats) # let's us easily calculate stats on rows or columns of a data matrix
library(cowplot) # allows you to combine multiple plots in one figure
library(gt)
library(RColorBrewer)
library(gplots)
#load study design file
targets <- read_tsv("data/studydesign.txt")

#Create paths to each sample
path <- file.path("results/mappedReads/",targets$sample, "abundance.tsv")
type(path)

#Check if the paths are valid
all(file.exists(path)) 

# get annotations using organism-specific package ----

#transcripts is a function of ensemble package

Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))

Tx <- as_tibble(Tx)

#need to change first column name to 'target_id'
Tx <- dplyr::rename(Tx, target_id = tx_id)

#transcript ID needs to be the first column in the dataframe
Tx <- dplyr::select(Tx, "target_id", "gene_name")

# import Kallisto transcript counts into R using Tximport ----

# copy the abundance files to the working directory and rename so that each sample has a unique name

Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx, #It should specifically have target_id first and then gene_name
                     txOut = FALSE, #gene level summary
                     countsFromAbundance = "lengthScaledTPM", 
                     ignoreTxVersion = TRUE) #ignores the.1,.2 versions of the transcript ids

#If you want the number of transcripts
Txi_transcripts <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx, #It should specifically have target_id first and then gene_name
                     txOut = TRUE, #It gives transcript level summary
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE) #ignores the.1,.2 versions of the transcript ids


#exploring the data to find if the total counts in transcripts,gene and counts are making sense
myTPM <- Txi_gene$abundance
myTPM_transcrits <- Txi_transcripts$abundance
myCounts <- Txi_gene$counts


#filling in the sample ids so that they can be given to column names easily
sampleLabels <- targets$sample

# Apply Filtering and Normalization----
#DGE list is a function of edgR package that converts the matrix into a list which contains normalization factors,count matrix and sample metadata
myDGEList <- DGEList(Txi_gene$counts)

#log 2 normalization is done to remove heteroscadacity and for visualizing our data
log2.cpm <- cpm(myDGEList, log=TRUE)

log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
colnames(log2.cpm.df) <- c("geneID", sampleLabels)

#making the data tidy for visualization
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, 
                                  cols = HS01:CL13, 
                                  names_to = "samples", 
                                  values_to = "expression") 


# Plot to clearfy the effect of TMM normalization and filtering low cpm data 
p1 <- ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

#remove the outliers that are present in more than five samples and visualize
cpm <- cpm(myDGEList)
keepers <- rowSums(cpm>1)>=5
myDGEList.filtered <- myDGEList[keepers,]

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, 
                                           cols = HS01:CL13, 
                                           names_to = "samples", 
                                           values_to = "expression") 

p2 <- ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, # dataframe to be pivoted
                                                cols = HS01:CL13, # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)

p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

combined_plot <- plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)

ggsave("results/plots/CPM_violin_TMM_filtering_comparison.png", 
       plot = combined_plot,
       width = 12, height = 8, dpi = 300, bg= "white")


#Dimensionality reduction ----

# Plot PCA result to clarify two groups
group <- targets$group
group <- factor(group)


pca.res <- prcomp(t(log2.cpm.filtered.norm), scale. = FALSE, retx = TRUE)#scale is false because we already did normalization and transformation
summary(pca.res)
pc.var <- pca.res$sdev^2
pc.per <- round(pc.var / sum(pc.var) * 100, 1)


pca.res.df <- as_tibble(pca.res$x)
pca.res.df$sra_accession <- targets$sra_accession
pca.res.df$sample <- targets$sample
pca.res.df$group <- factor(targets$group)

pca.plot_sra_accession <- ggplot(pca.res.df) +
  aes(x = PC1, y = PC2, color = group, label = sra_accession) +
  geom_point(size = 4) +
  geom_text_repel(size = 3.5, show.legend = FALSE) +
  stat_ellipse() +
  xlab(paste0("PC1 (", pc.per[1], "%)")) +
  ylab(paste0("PC2 (", pc.per[2], "%)")) +
  labs(title = "PCA Plot: TMM-Normalized Expression",
       caption = paste0("Produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

pca.plot_sample <- ggplot(pca.res.df) +
  aes(x = PC1, y = PC2, color = group, label = sample) +
  geom_point(size = 4) +
  geom_text_repel(size = 3.5, show.legend = FALSE) +
  stat_ellipse() +
  xlab(paste0("PC1 (", pc.per[1], "%)")) +
  ylab(paste0("PC2 (", pc.per[2], "%)")) +
  labs(title = "PCA Plot: TMM-Normalized Expression",
       caption = paste0("Produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

ggsave("results/plots/PCA_TMM_disease_VS_healthy.png", 
       plot = pca.plot_sample,
       width = 12, height = 8, dpi = 300, bg= "white")



#differential expression analysis----

# setting up the study design - linear model with our target variables

group <- factor(targets$group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

#model the mean variance relationship
#voom function of limma package will model the mean variance relationship

v.DGEList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = TRUE)

#fit a normal model to your data
fit <- lmFit(v.DGEList.filtered.norm,design)

# Contrast matrix 
contrast.matrix <- makeContrasts(infection = disease - healthy,
                                 levels=design)

# extract the linear model fit 
fits <- contrasts.fit(fit, contrast.matrix)
#get bayesian stats for your linear model fit
ebFit <- eBayes(fits)

# TopTable to view DEGs 
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")

# convert to a tibble
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")


# Plot Venn diagram of the DEGs
vplot <- ggplot(myTopHits.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "Cutaneous leishmaniasis",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

ggplotly(vplot)

ggsave("results/plots/Vene_diagram.png", 
       plot = vplot,
       width = 12, height = 8, dpi = 300, bg= "white")

results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=1)
colnames(v.DGEList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DGEList.filtered.norm$E[results[,1] !=0,]
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")

# Cluster genes and samples for the heatmap
myheatcolors <- rev(brewer.pal(name="RdBu", n=11))
clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") #cluster rows by pearson correlation
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 

png("results/plots/heatmap_all_genes.png", width = 12, height = 8, units = "in", res = 300)
heatmap.2(diffGenes, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=myheatcolors, scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20))

dev.off()

# heatmap for second cluster (upregulated genes) 
modulePick <- 2 
myModule_up <- diffGenes[names(module.assign[module.assign %in% modulePick]),] 
hrsub_up <- hclust(as.dist(1-cor(t(myModule_up), method="pearson")), method="complete") 

png("results/plots/heatmap_upregulated_genes.png", width = 12, height = 8, units = "in", res = 300)
heatmap.2(myModule_up, 
          Rowv=as.dendrogram(hrsub_up), 
          Colv=NA, 
          labRow = NA,
          col=myheatcolors, scale="row", 
          density.info="none", trace="none", 
          RowSideColors=module.color[module.assign%in%modulePick], margins=c(8,20))
dev.off()


# heatmap for first cluster (downregulated genes) 
modulePick <- 1 
myModule_down <- diffGenes[names(module.assign[module.assign %in% modulePick]),] 
hrsub_down <- hclust(as.dist(1-cor(t(myModule_down), method="pearson")), method="complete") 

png("results/plots/heatmap_downregulated_genes.png", width = 12, height = 8, units = "in", res = 300)
heatmap.2(myModule_down, 
          Rowv=as.dendrogram(hrsub_down), 
          Colv=NA, 
          labRow = NA,
          dendrogram = "row",
          col=myheatcolors, scale="row", 
          density.info="none", trace="none", 
          RowSideColors=module.color[module.assign%in%modulePick], margins=c(8,20))
dev.off()

