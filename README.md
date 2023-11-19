# major_project

## Overview
The major project asseses a subset of RNA-Seq data obtained from a gene expression study https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182617.
The subset of data being analysed includes Col 0 ecotype _Arabidopsis thaliana_ plants that were trated with Salicylic acid at timepoint 0h and 24h

Following replicates are present 
File name | Treatment hour | FIle description 
--- | --- | ------- |
Col0_SA_treated_0h_rep1.fastq.gz | 0h | Biological replicate 1 (rep1) of Arabidopsis Col-0 ecotype (Col0) treated with Salicylic Acid (SA) for 0 hours (0h)
Col0_SA_treated_0h_rep2.fastq.gz | 0h | Biological replicate 2 (rep2) of Arabidopsis Col-0 ecotype (Col0) treated with Salicylic Acid (SA) for 0 hours (0h)
Col0_SA_treated_0h_rep3.fastq.gz | 0h | Biological replicate 3 (rep3) of Arabidopsis Col-0 ecotype (Col0) treated with Salicylic Acid (SA) for 0 hours (0h)
Col0_SA_treated_24h_rep1.fastq.gz | 24h | Biological replicate 1 (rep1) of Arabidopsis Col-0 ecotype (Col0) treated with Salicylic Acid (SA) for 24 hours (24h)
Col0_SA_treated_24h_rep2.fastq.gz | 24h | Biological replicate 2 (rep2) of Arabidopsis Col-0 ecotype (Col0) treated with Salicylic Acid (SA) for 24 hours (24h)
Col0_SA_treated_24h_rep3.fastq.gz | 24h | Biological replicate 3 (rep3) of Arabidopsis Col-0 ecotype (Col0) treated with Salicylic Acid (SA) for 24 hours (24h)

The final research goal was to identify the differentially expressed genes within the Salicylic acid hormone pathwayand different pathways regulated by salicylic acid.


---
title: "DEanalysis"
author: "Chano Ranwala"
date: "2023-11-16"
output: html_document
---

# Set up of libraries

As shown below the libraries required for this was set up, including providing a black an white global plot setting. Note that throughout the Rmarkdown `echo = FALSE` parameter was not added to the code chunk to as I wanted the printing of the R code that helped generate the results.


```{r warning=FALSE, message=FALSE}
# packages for DE analysis
library(tidyverse)
library(limma)
library(edgeR)

# packages for plot
library(ggplot2)
library(ggrepel)
library(reshape2)
library(scales)
library(ggsci)
library(Glimma)

# clear workspace
rm(list = ls())

# global plot setting
theme_set(theme_bw())
```

# Getting the gene annotation information

Gene annotation information was uploaded as shown below. This code below made the gene_anno_df data table that allows us to get the gene annotation information. Duplicates have been removed and the row names have been set. Model_name correlates to gene name. It provides us with a decription of the genes as well.


```{r}
#Get gene annotation information
gene_anno_df <- read.delim("https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_functional_descriptions")
gene_anno_df$Model_name <- gsub("\\..+", "", gene_anno_df$Model_name)
gene_anno_df <- gene_anno_df[!duplicated(gene_anno_df$Model_name), ]
rownames(gene_anno_df) <- gene_anno_df$Model_name

```

As shown below, the star_dir file and count_files were made. We basically set the Star dircetory file patrh. We also retrieved all files  that end with ReadsPerGene.out.tab which means all teh different replicates were retrived. 

```{r}
star_dir <- file.path("~/major_project/04_results/03_aligned_data")
count_files <- list.files(path = star_dir, pattern = "ReadsPerGene.out.tab$", full.names = T)
```

A raw count of all the reads for each gene all the different replicates at 0h and 24 h were made. The final raw_count_mt data matrix displays the gene name (column 1) as the unique identified for each preceeding column. Each preceeding column ( columns 2-6) shows the reeds observed per gene per condition.

```{r}
firstFile_df <- read.delim(count_files[1], header = F)
firstFile_df <- firstFile_df[-c(1:4), ] 
raw_count_mt <- matrix(0, nrow(firstFile_df), length(count_files))
rownames(raw_count_mt) <- firstFile_df[, 1]
colnames(raw_count_mt) <- rep("temp", length(count_files))

#This made the firstFile_df and raw_count_mt
for(i in 1:length(count_files)){
  sample_name <- gsub(".+\\/", "", count_files[i])
  sample_name <- gsub("\\..+", "", sample_name)
  tmp_df <- read.delim(count_files[i], header = F)
  raw_count_mt[,i] <- tmp_df[-c(1:4), 2]
  colnames(raw_count_mt)[i] <- sample_name
}
  
#This made the firstFile_df and raw_count_mt
count_genes <- gene_anno_df[rownames(raw_count_mt),] 
count_genes$gene_name <- rownames(raw_count_mt)
rownames(count_genes) <- count_genes$gene_name

#This made the tmp_df file. brought in variable i in values
raw_count_df <- cbind(gene_name = rownames(raw_count_mt), raw_count_mt)
write.csv(raw_count_df, file = "~/major_project/04_results/04_DE/Table1_reference_gene_raw_count.csv", row.names = F)

#count_genes file was made
```



raw_count_df made alongside the Table1_reference_gene_raw_count.csv file in the 04_DE folder

```{r}
counts <- raw_count_mt
samples <- colnames(raw_count_mt)
groups <- gsub("_rep.+$", "", colnames(raw_count_mt))
genes <- count_genes
dgeList<- DGEList(counts = raw_count_mt, 
    samples = samples,
    group = groups, 
    genes = count_genes)
```

countrs made in Data. genes, made in data. Samokes and groups made in values

DGE list made!

```{r}
dgeList$samples %>%
  rownames_to_column("sample") %>%
  ggplot(aes(x = sample, y = lib.size / 1e6, fill = group)) +
  geom_bar(stat = "identity") +
  ylab("Library size (millions)") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```
this chart was made

```{r}
dgeList <- calcNormFactors(dgeList)
```


something happened here god knows what it was

```{r}
cols = pal_npg("nrc", alpha=1)(2)
plotMDS(dgeList, labels = dgeList$samples$type, col = cols[dgeList$samples$group])
```

```{r}
head(cpm(dgeList))
```



```{r}
plotDensities(cpm(dgeList, log = TRUE), main = "all genes", legend = "topright")
genes2keep <- rowSums(cpm(dgeList) > 1) > 3
summary(genes2keep)
```

chart madeeee

```{r}
plotDensities(cpm(dgeList, log = TRUE)[genes2keep, ], main = "fitered genes", legend = "topright")
```



```{r}
dgeList <- dgeList[genes2keep, ]
```


```{r}
design <- model.matrix(~0 + group, data = dgeList$samples)
design
dgeList <- estimateDisp(dgeList, design = design)
```


```{r}
etResults <- exactTest(dgeList, pair = 1:2)$table
hist(etResults[,"PValue"], breaks= 50)
```


```{r}
etResults <- etResults %>%
  rownames_to_column("Geneid") %>%
  mutate(FDR = p.adjust(PValue)) %>%
  arrange(PValue) %>%
  as_tibble()
```

```{r}
positive_count <- sum(etResults[, 2] > 0)
negative_count <- sum(etResults[, 2] < 0)
# Print the counts
cat("Positive count:", positive_count, "\n")
cat("Negative count:", negative_count, "\n")
```


This was for me to count how many genes were positively expressed and which were negatively expressed




```{r}
cpm(dgeList, log= TRUE)["AT5G48720",]
```

```{r}
sigGenes <- filter(etResults, FDR< 0.01, abs(logFC) > 1)$Geneid
```

```{r}
write.csv(etResults[etResults$Geneid %in% sigGenes, ], 
            file = "~/major_project/04_results/04_DE/Table2_sigGenes.csv", 
            row.names = F)
result_df <- read.csv("~/major_project/04_results/04_DE/Table2_sigGenes.csv", header = TRUE)

# Count the number of significant genes (excluding the header row)
num_sig_genes <- nrow(result_df)

# Print the number of significant genes
cat("Number of significant genes:", num_sig_genes)
```

```{r}
plotMD(dgeList, status = rownames(dgeList) %in% sigGenes)
```

```{r}
etResults %>%
  mutate(DE = Geneid %in% sigGenes) %>%
  ggplot(aes(x = logFC, y = -log10(PValue),colour = DE)) +
  geom_point() +
  scale_colour_manual(values = c("grey50", "red"))
```

```{r}
DEs <- decideTests(exactTest(dgeList), adjust.method = "fdr", p.value = 0.01, lfc = 1)
glMDPlot(exactTest(dgeList), counts = cpm(dgeList, log = TRUE),
        groups = dgeList$samples$group,
        status = DEs[, 1],
        main="MD plot: Col_0_normal_light_vs_Col_0_mock",
        side.main = "ID", side.ylab = "Expression (CPM)",
        sample.cols = pal_npg()(2)[dgeList$samples$group],
        folder = "glimmer_plot", launch = FALSE)
```

```{r}

```
