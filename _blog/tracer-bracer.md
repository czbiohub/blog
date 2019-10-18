---
layout: post
mathjax: true
title: Tools to Explore the Immune Repertoire
date: 2019-10-20
author: Clarissa Vazquez-Ramos
---

We are interested in how the immune system changes with age, and we accomplish that by building a clonal landscape of T and B cells to illustrate the diversity of the immune repertoire.


## What are T and B Cells?
<img align="right" src="/images/tracer-bracer/bcell.png" width="25%" height="25%">

T and B cells are lymphocytes and active participants of the adaptive immune response. These cells use the molecules on their surface, T cell receptors (TCR) and B cell receptors (BCR), to recognize pathogens in the body, such as a virus. T and B cells undergo V(D)J recombination — a process in which their DNA is shuffled — in order to develop their receptor which can potentially recognize pathogens. TCRs and BCRs accomplish this with the unique biological code they express in which can identify only a *specific antigen*; it can be thought of as two puzzle pieces meant to fit together. Although T and B cells serve the same purpose as participants in the adaptive immune response, the difference therein lies in their interactions with antigen-presenting cells and clonal expansion.

A T cell attacks and kills an antigen-presenting cell directly. When the B cell receptor recognizes its cognate antigen, it secretes several free forms of its receptor (antibodies). These antibodies are what latch onto the antigen and mark it for destruction by other cells.

<p align="center">
<img src="/images/tracer-bracer/tcell-bcell.png" alt="tcell vs bcell" width="45%" height="45%">
</p>

When T and B cells encounter antigen-presenting cells, they undergo clonal expansion — a process in which T and B cells form clones of themselves to fight off pathogens, and remember them for future encounters. It is important to note though that when B cells proliferate, their receptors undergo hypermutation — a process in which BCRs are diversified to enable the immune system to adapt its responses to new threats. These mutations occur at an extremely high rate, so during clonal expansion thousands of B cells may possess slightly different receptors. The B cell with the highest affinity towards the antigen will be selected to produce antibodies.

We can explore clonal expansion in single-cell datasets with the tools TraCeR and BraCeR.

## TraCeR & BraCeR Core Functions
<img align="right" src="/images/tracer-bracer/tracer-bracer-fn.png" alt="tracer-bracer core functions" width="30%" height="30%">

[TraCeR](https://github.com/teichlab/tracer) and [BraCeR](https://github.com/teichlab/bracer) were developed specifically to handle single-cell data. The tools serve to reconstruct the sequences of TCR and BCR genes and identify cells that have the same receptor sequence. The two modes that perform this are *assemble* and *summarize*.

1. ***Assemble*** is almost identical in both TraCeR and BraCeR. They both take paired-end scRNA-seq reads and reconstruct the TCR and BCR sequences. The reconstructed sequences are used to identify cells that have the same receptor sequence. Reconstruction is accomplished with the following steps: [alignment](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [de novo assembly](https://github.com/trinityrnaseq/trinityrnaseq/wiki), [IgBlast](http://www.ncbi.nlm.nih.gov/igblast/faq.html#standalone), and TCR and BCR [expression quantification](http://pachterlab.github.io/kallisto/). BraCeR takes an extra step to perform a [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) search before IgBlast.
For each cell, an output directory is created with output from Bowtie2, Trinity, (BLAST), IgBlast, and Kallisto as well as files describing the TCR and BCR sequences that were assembled.

2. ***Summarize*** takes the directories output from the ***assemble*** phase of several cells. Recall that B cells undergo hypermutation during clonal expansion, so the results in this step differ for TraCeR and BraCeR. In TraCeR, clone assignments are generated for each cell. In BraCeR, we receive a database file containing all the reconstructed sequences (e.g. CDR3, V and J).

For TraCeR, we don't need to take any further steps after `summarize`, as we have the clonal groups already assigned. For BraCeR, extra steps outside of the tool are necessary to generate clone assignments. Clone assignment is accomplished by first dividing the antibody heavy chain variable region (VH) sequences into groups which contain the same V and J genes and CDR3 length. A clone assignment is made if the amino acid CDR3 sequence shares *similarity* with other members within its groups.

Overall, these tools helps us identify single cells that have undergone clonal expansion. We currently run TraCeR and BraCeR analysis pipelines on [AWS Batch](https://docs.aws.amazon.com/batch/latest/userguide/what-is-batch.html) by manually submitting the jobs. We submit thousands of cells to be assembled asynchronously, then pull the assembled cells down and summarize them to identify clonal groups. While this workflow already carries out the analysis, we wanted to improve its reproducibility. Thus, we turned to Nextflow.


## What is Nextflow?
[Nextflow](https://www.nextflow.io/) is a workflow manager which allows for scalable and reproducible scientific workflows using containers. It simplifies the implementation and deployment of complex, parallel workflows, which is necessary for the thousands of single cells we process. Because Nextflow is based on the dataflow programming model, you can effortlessly link processes together in one workflow.

Nextflow has the capability to run pipelines in the cloud (e.g. AWS Batch) or locally.


## The Implementation
The implementation performs three steps which are linked together through `Channels`. A `Channel` has two major properties: sending messages and receiving data. A `Channel` sends messages in an asynchronous manner in which the operation will complete immediately, without having to wait for the receiving process. It will also receive data, which is a blocking operation where the receiving process is stopped until the message has arrived. The figure below illustrates how the workflow manages the three different steps. Steps 1 and 2 are handled asynchronously while Step 3 relies on the completion of Step 2 before starting.
<p align="center">
<img src="/images/tracer-bracer/nf-workflow.png" alt="nextflow workflow" width="70%" height="70%">
</p>


### Step 1: Preparation
The first step is to prepare the reads for the next processes. In this specific workflow, zipped fastq pair files are expected. Unzipping them first is necessary as the next steps only work wit unzipped files.

**`process unzip_reads:`**
> The first step is to open a `Channel` using the method [`.fromFilePairs()`](https://www.nextflow.io/docs/latest/channel.html#fromfilepairs). This method returns the file pairs matching the [glob](https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob) pattern input by the user. The fastq files get unzipped and passed into a new `Channel`, ***reads_unzipped_ch***,for the following process.


### Step 2: Assembly
In this step, we assemble the reads with either TraCeR or BraCeR. The figure below illustrates a run with TraCeR.

**`process assemble:`**
> This process takes in the unzipped fastq files from `reads_unzipped_ch` and reconstructs the TCR or BCR sequences, depending which analysis is being ran. The reads are assembled asynchronously and the output folders are published to a user-specified directory. These same folders are also passed into a new `Channel`, ***assembled_ch***, which will be used for the last process.


### Step 3: Summarize
Finally, in this last step we generate summary statistics and begin clone assignment.

**`process summarize:`**
> This last process calls the method [`.collect()`](https://www.nextflow.io/docs/latest/operator.html#operator-collect) on ***assembled_ch***. This method collects all the files emitted from ***assembled_ch*** into a list and uses that as input for this process. The output is published to the same directory where the assembled files are.

![nf-tracer run](/images/tracer-bracer/nf-tracer.gif)


### Step 4: Visualization
Different visualizations you could create with output from TraCeR/BraCeR are clonal trees, clonal networks, pie charts, etc. A clonal network can help us visualize the landscape a cell population's clonal groups. For example, in the figures below we visualize the landscape of our [Tabula Muris Senis](https://github.com/czbiohub/tabula-muris-senis) dataset by observing the number of clones for two age groups: 3 months and 24 months. This network shows us that the number of clones has increased with age. If we look at the distribution of clonal cells vs. singletons (non-clonal cells) through a pie chart instead, we see that the ratios of singletons to clonal cells have changed.

<img align="right" src="/images/tracer-bracer/bracer-output.png" alt="bracer output" width="48%" height="48%">
<img src="/images/tracer-bracer/tracer-output.png" alt="tracer output" width="48%" height="48%">
<p align="center">
<img src="/images/tracer-bracer/legend.png" alt="legend" width="50%" height="50%">
</p>
