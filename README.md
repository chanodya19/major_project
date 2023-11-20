# major_project

# Overview
The major project asseses a subset of RNA-Seq data obtained from a gene expression study https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182617.
The subset of data being analysed includes Col 0 ecotype _Arabidopsis thaliana_ plants that were treated with Salicylic acid at timepoint 0h and 24h post treatment.

Following replicates are present. 
File name | Treatment hour | File description 
--- | --- | ------- |
Col0_SA_treated_0h_rep1.fastq.gz | 0h | Biological replicate 1 (rep1) of Arabidopsis Col-0 ecotype (Col0) treated with Salicylic Acid (SA) for 0 hours (0h)
Col0_SA_treated_0h_rep2.fastq.gz | 0h | Biological replicate 2 (rep2) of Arabidopsis Col-0 ecotype (Col0) treated with Salicylic Acid (SA) for 0 hours (0h)
Col0_SA_treated_0h_rep3.fastq.gz | 0h | Biological replicate 3 (rep3) of Arabidopsis Col-0 ecotype (Col0) treated with Salicylic Acid (SA) for 0 hours (0h)
Col0_SA_treated_24h_rep1.fastq.gz | 24h | Biological replicate 1 (rep1) of Arabidopsis Col-0 ecotype (Col0) treated with Salicylic Acid (SA) for 24 hours (24h)
Col0_SA_treated_24h_rep2.fastq.gz | 24h | Biological replicate 2 (rep2) of Arabidopsis Col-0 ecotype (Col0) treated with Salicylic Acid (SA) for 24 hours (24h)
Col0_SA_treated_24h_rep3.fastq.gz | 24h | Biological replicate 3 (rep3) of Arabidopsis Col-0 ecotype (Col0) treated with Salicylic Acid (SA) for 24 hours (24h)

The final research goal was to identify the differentially expressed genes within the salicylic acid hormone pathway and different pathways regulated by salicylic acid.

# Differential Expression Analysis.
---
title: "DEanalysis"
author: "Chano Ranwala"
date: "2023-11-16"
output: html_document
---

## Set up of libraries

As shown below the libraries required for this was set up, including providing a black and white global plot setting. Note that throughout the Rmarkdown `echo = FALSE` parameter was not added to the code chunk to as I wanted the printing of the R code that helped generate the results.


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

## Getting the gene annotation information

Gene annotation information was uploaded as shown below. This code below made the gene_anno_df data table that allows us to get the gene annotation information. Duplicates have been removed and the row names have been set. Model_name correlates to gene name. It provides us with a description of the genes as well.


```{r}
#Get gene annotation information
gene_anno_df <- read.delim("https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_functional_descriptions")
gene_anno_df$Model_name <- gsub("\\..+", "", gene_anno_df$Model_name)
gene_anno_df <- gene_anno_df[!duplicated(gene_anno_df$Model_name), ]
rownames(gene_anno_df) <- gene_anno_df$Model_name

```

As shown below, the star_dir file and count_files were made. We basically set the Star directory file path. We also retrieved all files  that end with ReadsPerGene.out.tab which means all the different replicates were retrived. 

```{r}
star_dir <- file.path("~/major_project/04_results/03_aligned_data")
count_files <- list.files(path = star_dir, pattern = "ReadsPerGene.out.tab$", full.names = T)
```

A raw count of all the reads for each gene all the different replicates at 0h and 24 h were made. The final raw_count_mt data matrix displays the gene name (column 1) as the unique gene identified for each preceeding column. Each preceeding column ( columns 2-6) shows the reeds observed per gene per condition. 
raw_count_df table was made alongside the Table1_reference_gene_raw_count.csv file in the 04_DE folder

```{r}
#This made the firstFile_df and raw_count_mt
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
  
count_genes <- gene_anno_df[rownames(raw_count_mt),] 
count_genes$gene_name <- rownames(raw_count_mt)
rownames(count_genes) <- count_genes$gene_name

raw_count_df <- cbind(gene_name = rownames(raw_count_mt), raw_count_mt)
write.csv(raw_count_df, file = "~/major_project/04_results/04_DE/Table1_reference_gene_raw_count.csv", row.names = F)

```

Different variables were allocated to the global environment. DGE list construction was initialised, starting with the library sizes being shown in the form of number of reads per replicate.

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

```{r}
dgeList$samples %>%
  rownames_to_column("sample") %>%
  ggplot(aes(x = sample, y = lib.size / 1e6, fill = group)) +
  geom_bar(stat = "identity") +
  ylab("Library size (millions)") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

When sequencing the reads, we recieved different number of genes mapped to the genome for each replicate. Therefore we used the following command to normalise read counts based on the library of each sample. 

```{r}
dgeList <- calcNormFactors(dgeList)
```


## Statistical analysis based on DGEList

Here below we used a multiple dimentional scaling plot to check the overall gene expression across different samples. This generated the MDS plot.

```{r}
cols = pal_npg("nrc", alpha=1)(2)
plotMDS(dgeList, labels = dgeList$samples$type, col = cols[dgeList$samples$group])
```

Our samples our grouped correctly as per the generated plot. Different conditions ( 0h versus 24h) were grouped together. No different conditions were mixed together. A clear pattern of separation was noticed. 

After the overall gene expression was determined, we utilised a counts per milllion factor to have a look at the normalised read conda value for each gene.

```{r}
head(cpm(dgeList))
```
This gives the first 6 genes and how they are expressed within samples.

## Applying further adjustments and filtering

Filtrations of low expression henees is very important as it will allow the statistical power of the study to increase. The expression intensity of all genes in all samples were plotted against their density using the following command. Removed genes that express less than 1 intensity level (CPM ). A summary of the genes that fulfilled this condition was also checked

```{r}
plotDensities(cpm(dgeList, log = TRUE), main = "all genes", legend = "topright")
genes2keep <- rowSums(cpm(dgeList) > 1) > 3
summary(genes2keep)
```

The final plot with the filtered genes were checked through the command below. 

```{r}
plotDensities(cpm(dgeList, log = TRUE)[genes2keep, ], main = "fitered genes", legend = "topright")
```
Dataset was filltered that allowed only the genes with good abundances remained in dgeList.
```{r}
dgeList <- dgeList[genes2keep, ]
```

### Dispersion calculation 

The following design shows a matrix that states which samples belong to which group. The dispersion was then estimated, which allowed to proceed with differential expression analysis.

```{r}
design <- model.matrix(~0 + group, data = dgeList$samples)
design
dgeList <- estimateDisp(dgeList, design = design)
```

This below conducts a fishers exact test on the two 0h and 24 h group. Next a histogram was generated to see the number of genes that show which p value in the range from 0-1.0.

```{r}
etResults <- exactTest(dgeList, pair = 1:2)$table
hist(etResults[,"PValue"], breaks= 50)
```
We then adjusted the p value using False Discovery Rate (FDR). 

```{r}
etResults <- etResults %>%
  rownames_to_column("Geneid") %>%
  mutate(FDR = p.adjust(PValue)) %>%
  arrange(PValue) %>%
  as_tibble()
```

In the final etResults table, I wanted to see which genes were overexpressed and which genes were down regulated using the code below. This was for me to count results many genes were positively expressed and which were negatively expressed, but these are from all the filtered genes, these counts are not a representation of the finalised significant results. 

```{r}
positive_count <- sum(etResults[, 2] > 0)
negative_count <- sum(etResults[, 2] < 0)
# Print the counts
cat("Positive count:", positive_count, "\n")
cat("Negative count:", negative_count, "\n")
```


Finally, we make the table with significantly expressed genes. The FDR was set to 0.01 while the logFC was set to more than 1.
I also checked how many significant genes are present and which are upregulated/downregulated.

```{r}
sigGenes <- filter(etResults, FDR< 0.01, abs(logFC) > 1)$Geneid
positive_count <- sum(sigGenes[, 2] > 0)
negative_count <- sum(sigGenes[, 2] < 0)
# Print the counts
cat("Positive count:", positive_count, "\n")
cat("Negative count:", negative_count, "\n")

#write CSV file with significantly differentionally expressed genes
write.csv(etResults[etResults$Geneid %in% sigGenes, ], 
            file = "~/major_project/04_results/04_DE/Table2_sigGenes.csv", 
            row.names = F)
result_df <- read.csv("~/major_project/04_results/04_DE/Table2_sigGenes.csv", header = TRUE)

# Count the number of significant genes (excluding the header row)
num_sig_genes <- nrow(result_df)

# Print the number of significant genes
cat("Number of significant genes:", num_sig_genes)
```
The final CSV file with significant gene are now made. This is needed for further KEGG pathway analysis using DAVID.

Next the pattern of logFC to expression level was demonstrated through a chart.
```{r}
plotMD(dgeList, status = rownames(dgeList) %in% sigGenes)
```
I also made a volcano plot which plotted logFC values against p values as an alternative means of analysis.

```{r}
etResults %>%
  mutate(DE = Geneid %in% sigGenes) %>%
  ggplot(aes(x = logFC, y = -log10(PValue),colour = DE)) +
  geom_point() +
  scale_colour_manual(values = c("grey50", "red"))
```
## Glimma plot

Finally a Glimma plot folder was made. This allowed an interactive HTML plot to be generated that shows all genes together alongside their information, relative normalised expression value and expression at 0h versus 24h.

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

End.
