---
layout: post
mathjax: true
title: Tabula Muris
---

The Chan Zuckerberg Biohub recently released Tabula Muris, a compendium of single cell transcriptome data from the mouse containing nearly 100,000 cells from 20 organs and tissues. The data allow for direct and controlled comparison of gene expression in cell types shared between tissues, such as immune cells from distinct anatomical locations. They also allow for a comparison of two distinct technical approaches:

* microfluidic droplet-based 3'-end counting, which provides a survey of thousands of cells per organ at relatively low coverage.
* FACS-based full length transcript analysis, which provides higher sensitivity and coverage.

We hope this rich collection of annotated cells will be a useful resource for:

* Defining gene expression in previously poorly-characterized cell populations.
* Validating findings in future targeted single-cell studies.
* Developing of methods for integrating datasets (eg between the FACS and droplet experiments), characterizing batch effects, and quantifying the variation of gene expression in many cell types between organs and animals.

You can access the data at a few levels of depth:

* A simple [website](http://tabula-muris.ds.czbiohub.org/) shows the distribution of gene expression in each tissue and cell population.
* Matrices of gene-cell counts and metadata are available as CSVs on Figshare (and as direct downloads)
* The fastq files for each library are available on the [Short Read Archive](https://www.ncbi.nlm.nih.gov/sra/?term=SRP131661).

For some organs, like the Bladder, Diaphragm, and Tongue epithelium, this represents the first large single-cell sequencing study to our knowledge. Tabula Muris is also sex-balanced, in [constrast](https://www.the-scientist.com/?articles.view/articleNo/48616/title/How-Much-Do-Sex-Differences-Matter-in-Mouse-Studies-/) to the majority of murine research. (This bias is so strong that some [studies](http://www.cell.com/action/showMethods?pii=S2405-4712%2816%2930265-4) decline even to state the sex of the mice studied.) It represents the first large single-cell sequencing study of the female murine liver, kidney, and skin to our knowledge.

<img align="right" src="/images/2018-8-10-Introduction/list_of_organs.jpg" width="300" height="300">

One of the challenges of analyzing single-cell sequencing data is the role of batch effects. Tissues from genetically identical mice processed by different labs using different dissociation, library preparation, and sequencing protocols may produce rather different profiles of gene expression. To provide a cross-organ dataset which is as standardized as possible, the mice were processed in a rather delicate ballet: representatives from over a dozen labs (each handling an organ or two) would gather as the mouse was sacrificed, extract their lab's organ, dissociate it, and bring it to a central location for sorting and sequencing.

![jpg](/images/2018-8-10-Introduction/crowd for organs.jpg)
*A crowd waits for organs*

The payoff for taking that care can be seen in a TSNE plot of all the FACS-sorted cells together, where cells from different tissues mix together in visible clusters. (In fact, these regions represent common cell types such as epithelial cells, endothelial cells, immune cells, etc.)

![](/images/2018-8-10-Introduction/facs_tsne_by_tissue.png)
*TSNE of all FACS-sorted from all organs.*

## Analysis

Cell type annotation in single-cell sequencing is, at present, a semi-supervised problem. It requires expert knowledge (what marker genes are thought to be uniquely expressed in a given cell type) and high-dimensional data analysis (extracting information from a collection 20,000-dimensional gene expression vectors). To empower the organ experts from each of the collaborating labs to analyze the data they collected and to make the analysis legible to the community at large, we elected to use a relatively simple pipeline as instantiated in the software package Seurat.

The pipeline begins with vectors of gene expression in each cell, normalizes the data, selects highly variable genes, reduces their dimension using PCA, then clusters cells based on a nearest-neighbors graph. By inspecting the expression of marker genes (and occasionally computing differentially expressed genes when the markers were not present), researchers were able to identify each cluster. In the case that the profile of gene expression indicated a mixture of populations, they would either tune parameters (like number of PCs ) or subset it and repeat the whole pipeline. To see a worked example, check out the [Organ Annotation Vignette](/files/Organ_Annotation_Vignette.pdf), which also describes all of the mathematical functions and parameter values used. The 33 notebooks used to analyze each tissue and experimental method can be found in [github](https://github.com/czbiohub/tabula-muris/tree/master/00_data_ingest/02_tissue_analysis_rmd).

Note that this method does not require that one resolve batch effects; it is enough that, within each batch, cells of different types separate. For example, in the limb muscle, there were two clusters of cells (largely from different mice) that were both recognizable by gene expression as satellite cells.

<img src="/images/2018-8-10-Introduction/muscle.png" alt="muscle tsne" width="1500px"/>

To facilitate the reuse of this data, we provide annotations in the controlled vocabulary of a [cell ontology](http://obofoundry.org/ontology/cl.html) from The OBO Foundry. Since definitions of cell types are quickly evolving and ontologies require time to catch up, we also included a free annotation field where arbitrary notes could be added. Those include features like subtypes (luminal progenitor cells in the mammary gland) or localization (periportal hepatocytes).

## Vignettes

In the coming weeks, we will have a series of posts highlighting ways to use this resource. We will talk about:

* The impact of cell type heterogeneity on regressing confounders
* Building predictors from annotated data (and using them to clean the annotations you already have)
* Patterns of gene expression in B cells found in different organs in the body
* Variability in 'differentially expressed genes' found in three datasets, our smartseq2 and 10X data and the microwell-seq from the Mouse Cell Atlas of Han et Al.
* etc...

We would also like to invite you to contribute to this series, especially to showcase methods that address major problems in defining cell type/state, aligning large datasets, and finding latent structures.

Stay tuned!

-Josh Batson and the TM Consortium
