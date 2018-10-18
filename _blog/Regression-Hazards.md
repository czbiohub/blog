---
layout: post
mathjax: true
codehide: true
title: Regression Hazards
date: 2018-10-20
author: Josh Batson
---


tl;dr "Technical confounders" are likely correlated with real biology. Regressing them out can badly distort the data.

When analyzing single-cell sequencing data, it is natural to want to remove the effect of technical confounders. Statistics such as the cellular detection rate (number of genes detected per cell), library size (number of reads or UMI per cell), or percent mitochondrial reads may all represent some kind of bias. If those statistics vary between cell types, however, regressing them out will pull the expression vectors for those cell types closer together.

We can see a few examples of this phenomenon in the Tabula Muris dataset.

```r
library(Seurat)
library(here)
library(tidyverse)

tm.droplet.matrix = readRDS(here("data", "TM_droplet_mat.rds"))

tm.droplet.metadata = read_csv(here("data", "TM_droplet_metadata.csv"))
row.names(tm.droplet.metadata) = tm.droplet.metadata %>% pull(cell)

tmd <- CreateSeuratObject(raw.data = tm.droplet.matrix, meta.data = tm.droplet.metadata, project = "TabulaMuris")

# Only keep annotated cells
annotated_cells = tm.droplet.metadata %>% filter(!is.na(cell_ontology_class)) %>% pull(cell)
tmd <- SubsetData(tmd, cells.use = annotated_cells, do.clean = TRUE)

tmd <- NormalizeData(tmd)
```

```PlainText
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
```

Consider the Bladder. Because the epithelial cells have significantly more UMI on average than the mesenchymal or endothelial cells, we are set up for Simpson's paradox: even if the expression of a gene is positively correlated with nUMI within each cell type, it may be negatively correlated if all cell types are considered together.

```r
FetchData(tmd, c('nUMI', 'tissue', 'free_annotation'))%>% filter(tissue == 'Bladder') %>%
  ggplot(aes(x = tissue, y = log(nUMI), fill = free_annotation)) + geom_boxplot() + coord_flip()
```

![png](/images/Regression-Hazards/regression-hazards_0.png)

For example, the mesenchymal marker _Car3_ is positively correlated with nUMI within each cell type and negatively correlated at the level of the ensemble.

```r
FetchData(tmd, c('nUMI', 'tissue', 'Car3', 'free_annotation')) %>% filter(tissue == 'Bladder') %>%
  ggplot(aes(log(nUMI), Car3, color = free_annotation)) + geom_jitter(width = 0.2, alpha = 0.3) + geom_smooth(method = 'lm') +
  geom_smooth(method = 'lm', color = 'black')
```

![png](/images/Regression-Hazards/regression-hazards_1.png)

A similar phenomenon can be seen in the tongue, where the number of genes detected per cell differs between cell types.

```r
FetchData(tmd, c('nGene', 'tissue', 'cell_ontology_class'))%>% filter(tissue == 'Tongue') %>%
  ggplot(aes(x = tissue, y = nGene, fill = cell_ontology_class)) + geom_boxplot() + coord_flip()
```

![png](/images/Regression-Hazards/regression-hazards_2.png)

Here, the proliferation marker _Top2a_ has quite different correlations with nGene globally and within each cell type.

```r
FetchData(tmd, c('nGene', 'tissue', 'Top2a', 'cell_ontology_class')) %>% filter(tissue == 'Tongue') %>%
  ggplot(aes(nGene, Top2a, color = cell_ontology_class)) + geom_jitter(width = 0.2, alpha = 0.3) + geom_smooth(method = 'lm') +
  geom_smooth(method = 'lm', color = 'black')
```

![png](/images/Regression-Hazards/regression-hazards_3.png)

In fact, this sort of paradox is possible in every organ we collected. Each organ has a cell type which is an ourlier for nGene, setting up a potential Simpson's paradox for each gene differentially expressed in that cell type.

```r
FetchData(tmd, c('nGene', 'tissue', 'cell_ontology_class'))%>%
  ggplot(aes(x = tissue, y = nGene, fill = cell_ontology_class)) + geom_boxplot() + coord_flip() +
  ggtitle('nGene variability between cell types within tissues') +
  theme(axis.text = element_text(size=40), text = element_text(size=40), legend.position="none")
```

![png](/images/Regression-Hazards/regression-hazards_4.png)

Any form of regularization (regressing confounders, batch effect correction) has the potential to introduce artifacts. This is why it's essential that, whatever regularization you used to find clusters or trajectories, you use raw gene counts for computing differential expression. Terms like nUMI, nGene, or batch can be used as covariates in a regression model for gene expression. By giving that final model access to the covariates, instead of trying to 'adjust' for them beforehand, you increase the likelihood that the DE genes you find come from real biological variation.
