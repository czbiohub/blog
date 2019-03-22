---
layout: post
mathjax: true
title: Towards a phylogeny of cell types
date: 2019-03-22
author: Olga Botvinnik
---


[TOC]: # "Table of Contents"






## Introduction

- What features contribute to speciation?
- Does introduction of a new cell type create a new species?
- cells --> tissues --> organs --> organ systems --> organisms
- Cells are a fundamental unit of life
- To study evolution of species, need a quantitative method of studying phenotype
  - single-cell RNA-seq is a quantitative phenotypic readout

## Compare cell types across species with protein $k$-mers


### Methods

#### Compute a compressed k-mer representation of translated protein k-mers from each sample

figure: protein k-mer min hashing

##### Input(s)

- One pair of fastq.gz files per sample/single cell **already supported by `sourmash`**
- One bam file per 10x genomics run **added to `sourmash`, but is poorly implemented (appends to a Python list ..)**
  - One bam file per inDrop run **not yet supported by `sourmash`**

##### Output

One k-mer signature `.sig` file containing a sketch of protein k-mer hashes

#### Remove species-specific signal with TF-IDF


##### Input(s)

A protein $k$-mer signature file per sample and a csv matching the sample name
to the species label as a groupby method for removing the species signal.


##### Output

Species-removed k-mer sketch `.sig` files

#### Cluster samples on species-depleted protein k-mer similarity

1. Use Leiden or similar clustering
2. Potentially use "PAGA" or other single-cell style clustering


##### Input(s)

Species-removed protein $k$-mer signatures

##### Output(s)

1. Symmetric matrix of sample-sample $k$-mer similarities
2. A cluster label per cell to identify its cluster assignment
3. 

### How do we know we are right? "Gold standard" datasets

#### Recapitulate tissue-specificity from Merkin 2012 dataset

1. What about other datasets that cite this?

#### Hematopoeisis:

4. zebrafish, mouse, human, single-cell hematopoeisis single-cell datasets
5. Match Rosenthal 2018 tunicate dataset with mouse hematopoeisis, e.g. with Gene Expression Commons https://gexc.riken.jp/
6. Match enriched protein k-mers to bottom-up proteomics data for these cell types


### Questions we can ask with these methods


1. What peptide fragments are enriched in the different cell types? What genes do they match to?
2. What protein features are enriched in the different cells' k-mers? hydrophobicity, protein disorder, charge



### Software deliverables

Primarily I plan to use
[`sourmash`](https://github.com/dib-lab/sourmash/) to compute k-mer
similarities between samples

#### Potentially useful software rabbit holes

- Re-implement dashing/hyperloglog for both DNA and protein k-mers in Rust for `sourmash`
- Re-implement DOCKS/universal hitting k-mer set and only keep track of 21-mers
  fitting in length 80 seqs

#### Speed up I/O of signature saving/loading

Luiz Irber has
[started](https://github.com/dib-lab/sourmash/pull/532) this process by
[rewriting the backend](https://github.com/dib-lab/sourmash/pull/424) codebase of `sourmash` to [Rust](https://www.rust-lang.org/).

#### Add option for 3-frame translation (stranded RNA-seq data) to `sourmash`

#### Add option to perform TF-IDF to `sourmash`

1. e.g. as a `.csv` file containing a sample ID to cell type group mapping in
   `sourmash compute`
2. Or as a separate function `sourmash tfidf` on a set of sketches for storage
   and future use


- Methods to remove species-specific signal
  - Add option to perform within-group Term-frequency inverse document frequency to `sourmash`

#### Which clustering algorithms can be used to find groups of similar cell
types across species for a symmetric matrix?

  - Use symmetric similarity matrix as input to `scanpy`

#### Visualization of cell type clusters

  - PAGA, UMAP, nearest neighbors graph


## Build a phylogenetic tree of cell types, starting with hematopoeisis


### Methods

1. Download RNA-seq data across all species for hematopoeisis from SRA, ERA
2. Create groups of similar cells using **protein** k-mer simlarity
   1. Cross-check with published annotations from each paper
3. Calculate **nucleotide**-level k-mer similarities across all samples
4. Fix the cell type groups, and use nucleotide-level similarities across samples to build a tree using [Felstein's tree-pruning algorithm](https://en.wikipedia.org/wiki/Felsenstein%27s_tree-pruning_algorithm)

### How do we know we are right?

1. Use a species-specific cell type or tissue and ensure its "age" corresponds to the correct "age" of the species? - Questions

### Questions we can ask using this method

  1. Do "newer" cell types have more intrinsically disordered protein
     regions?

### Software deliverables

#### Fix memory management of 10x genomics input to `sourmash compute`

Currently appends all reads to Python list and runs out of memory. Should sort by cell barcode or similar to make only one barcode's signature at a time and avoid looking at ALL the data at once

### Add support of inDrop data input to `sourmash compute`

I've never worked with inDrop data before but there's a lot of single-cell data out there in this format

