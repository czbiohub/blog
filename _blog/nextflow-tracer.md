---
layout: post
mathjax: true
title: TraCeR & BraCeR
date: 2019-10-19
author: Clarissa Vazquez-Ramos
---

We are interested in how the immune system changes with age, and we accomplish that by building a clonal landscape of T and B cells to illustrate the diversity of the immune repertoire.


## What are T and B Cells?
<img align="right" src="/images/tracer-bracer/bcell.png" width="25%" height="25%">
T and B cells are lymphocytes (white blood cells) and active participants of the adaptive immune response. These cells serve to recognize pathogens through molecules found on their surface: T cell receptors (TCR) and B cell receptors (BCR). T and B cells undergo V(D)J recombination—a process in which their DNA is shuffled—to develop receptors and potentially recognize foreign entities in the body, such as viruses. TCRs and BCRs expresses unique biological code that recognize specific antigens. It can be thought of as a puzzle, where only one side of a puzzle piece is designed to fit with another. The difference between T and B cells lies in how they interact with external threats and antigen-presenting cells.

A T cell will attack and kill an antigen presenting cell directly. A B cell will secrete antibodies, leading it to turn into a plasma cell, and eventually lysing and releasing antibodies into the bloodstream.

<p align="center">
<img src="/images/tracer-bracer/tcell-bcell.png" alt="tcell vs bcell" width="40%" height="40%">
</p>

When T and B cells encounter antigen presenting cells, they undergo clonal expansion. Clonal expansion is a process in which T and B cells will multiply in order to fight off antigen presenting cells. We are able explore clonal expansion on our cells by using the tools TraCeR & BraCeR.

## TraCeR & BraCeR Core Functions
<img align="right" src="/images/tracer-bracer/tracer-bracer-fn.png" alt="tracer-bracer core functions" width="28%" height="28%">

[TraCeR](https://github.com/teichlab/tracer) and [BraCeR](https://github.com/teichlab/bracer) were developed specifically to handle single cell data. Their main purpose is to reconstruct the sequences of TCR and BCR genes and identify cells that have the same receptor sequence. The two modes that perform this are *assemble* and *summarize*

1. ***Assemble*** is almost identical in both TraCeR and BraCeR. They both take paired-end scRNA-seq reads and reconstruct their TCR and BCR sequences. The reconstructed sequences are used to identify cells that have the same receptor sequence. Reconstruction is accomplished with the following steps: [alignment](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [de novo assembly](https://github.com/trinityrnaseq/trinityrnaseq/wiki), [IgBlast](http://www.ncbi.nlm.nih.gov/igblast/faq.html#standalone), and TCR and BCR [expression quantification](http://pachterlab.github.io/kallisto/). BraCeR takes an extra step to perform a [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) search before IgBlast.
For each cell, an output directory is created with output from Bowtie2, Trinity, (BLAST), IgBlast, and Kallisto as well as files describing the TCR and BCR sequences that were assembled.

2. ***Summarize*** takes the directories output from the ***assemble*** phase of several cells, summarizes the receptor recovery rates, and generates clonotype networks from the assembled reads. This step helps us identify cells that have undergone clonal expansion.

Currently, we run TraCeR and BraCeR on [AWS Batch](https://docs.aws.amazon.com/batch/latest/userguide/what-is-batch.html) by manually submitting the jobs. We submit thousands of cells to be assembled asynchronously, then we pull the assembled cells down and summarize them to identify clonal groups. While this method already accomplishes our goal, we wanted to find a way to improve the reproducibility of the pipeline. Thus, we found Nextflow.


## What is Nextflow?
[Nextflow](https://www.nextflow.io/) is a workflow manager which allows for scalable and reproducible scientific workflows using containers. It simplifies the implementation and deployment of complex, parallel workflows, which is necessary for the thousands of single cells we process. Because Nextflow is based on the dataflow programming model, you can effortlessly link processes together in one workflow.

Nextflow has the capability to run pipelines on AWS Batch, but pipelines can also be ran locally.


## The Implementation
The implementation performs 3 steps which are linked together through `Channels`. A `Channel` has 2 major properties: sending messages and receiving data. A `Channel` sends messages in an asynchronous manner in which the operation will complete immediately, without having to wait for the receiving process. It will also receive data, which is a blocking operation where the receiving process is stopped until the message has arrived. The figure below illustrates how the workflow manages the 3 different steps. Steps 1 and 2 are handled asynchronously while Step 3 relies on the completion of Step 2 before starting.
<p align="center">
<img src="/images/tracer-bracer/nf-workflow.png" alt="nextflow workflow" width="70%" height="70%">
</p>


### Step 1: Preparation
The first step is to open a `Channel` using the method [`.fromFilePairs()`](https://www.nextflow.io/docs/latest/channel.html#fromfilepairs). This method returns the file pairs matching the [glob](https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob) pattern input by the user. In our case, the file pairs returned are the sample names and the fastq files. These file pairs are then used as input for the first process.

**`process unzip_reads:`**
> In this first process, we prepare our fastq files for the next steps. We take the fastq files from the `Channel` we opened and unzip them. The unzipped fastq files and their respective sample names are passed into a new `Channel`, ***reads_unzipped_ch***, as output to be used in the following process.


### Step 2: TraCeR/BraCeR Assembly
In this step, we assemble the reads.

**`process assemble:`**
> This process takes in the unzipped fastq files from `reads_unzipped_ch` and reconstructs the TCR/BCR sequences. The reads are assembled asynchronously and the output is published to a specified directory. The directory is then output into a new `Channel`, ***assembled_ch***, which will be used for the last step.


### Step 3: TraCeR/BraCeR Summarize
Finally, in this last step we summarize the TCR and BCR recovery rates as well as generate clone assignments for the TCRs.

**`process summarize:`**
> This last process calls the method [`.collect()`](https://www.nextflow.io/docs/latest/operator.html#operator-collect) on ***assembled_ch***. What this does is collects all the files emitted from ***assembled_ch*** into a list and uses that as the input for `tracer/bracer summarize`. The output contains summary statistics describing successful TCR/BCR reconstruction rates as well information on the cells and which clonal groups they belong to. The output is published to the same directory where the assembled files are.

![nf-tracer run](/images/tracer-bracer/nf-tracer.gif)


### Step 4: Visualization
Different visualizations you could create with output from TraCeR/BraCeR are clonal trees, clonal networks, pie charts, etc. A clonal network can help us visualize the landscape a cell population's clonal groups. For example, in the figures below we visualize the landscape of our [Tabula Muris Senis](https://github.com/czbiohub/tabula-muris-senis) data by observing the number of clones for two age groups: 3 months and 24 months. This network shows us that the number of clones has increased with age. If we look at the distribution of clonal cells vs. singletons (non-clonal cells) through a pie chart instead, we see that the ratios of singletons to clonal cells have changed.

<img align="right" src="/images/tracer-bracer/bracer-output.png" alt="bracer output" width="48%" height="48%">
<img src="/images/tracer-bracer/tracer-output.png" alt="tracer output" width="48%" height="48%">
<p align="center">
<img src="/images/tracer-bracer/legend.png" alt="legend" width="50%" height="50%">
</p>
