---
layout: post
mathjax: true
title: TraCeR Pipeline: Integrated into a Nextflow workflow
date: 2019-09-06
author: Clarissa Vazquez-Ramos
---


## A Brief Introduction to TraCeR
[TraCeR](https://github.com/Teichlab/tracer) is a tool that reconstructs the sequences of rearranged and expressed T cell receptor genes from single-cell RNA-seq data. The TCR sequences are then used to identify cells that have the same receptor sequences and indicate that they are derived from the same original clonally-expanded cell.

Our TraCeR pipeline is currently ran from command line prompts and on AWS Batch. While this method already efficiently runs the pipeline, there is no straight forward approach to run other of the TraCeR tasks consecutively; Nextflow solves this problem.

## What is Nextflow?
[Nextflow](https://www.nextflow.io/) allows for scalable and reproducible scientific workflows using containers. It simplifies the implementation and deployment of complex, parallel workflows. Because Nextflow is based on the dataflow programming model, you can effortlessly link processes together in one workflow.

Nextflow also has the capability to run pipelines on AWS Batch without having to deal with the AWS interface.


## The Implementation
The implementation performs 3 tasks which are linked together through `Channels`. A `Channel` has 2 major properties: sending messages and receiving data. A `Channel` sends messages in an asynchronous manner in which the operation will complete immediately, without having to wait for the receiving process. It will also receive data, which is a blocking operation where the receiving process is stopped until the message has arrived.


### Step 1: Open a `Channel` and Preparation
The first step is to create a `Channel` using the method [`.fromFilePairs()`](https://www.nextflow.io/docs/latest/channel.html#fromfilepairs). This method returns the file pairs matching the [glob](https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob) pattern input by the user. In our case, the file pairs returned are the sample names and their fastq files of paired-end RNA-seq reads from single-cell. These file pairs are then used as input for the first process.

`process unzip_reads`:
In this first process, we prepare our fastq files for the next steps. We take the fastq files from the `Channel` we opened and unzip them. The unzipped fastq files and their respective sample names are passed into a new `Channel`, `reads_unzipped_ch`, as output to be used in the following process.


### Step 2: Assembly
In this step, we assemble the reads using TraCeR.

`process assemble`:
This process takes in the unzipped fastq files from `reads_unzipped_ch` and reconstructs the TCR sequences. The reads are assembled asynchronously and the output is published to a specified `S3 Bucket`. The bucket contains subdirectories for each sample with the output from Bowtie2, Trinity, IgBlast, Kallisto and Salmon as well as files describing the TCR sequences that were assembled.

The path to the directory containing all the above information is output into a new `Channel`, `assembled_ch`, which will be used for the last step.


### Step 3: Summarize
Finally, in this last step we summarize the TCR recovery rates as well as generate clonotype networks from the assembled reads.

`process summarize`:
This last process calls the method [`.collect()`](https://www.nextflow.io/docs/latest/operator.html#operator-collect) on `assembled_ch`. What this does is it collects all the files emitted from `assembled_ch` into a list and uses that as the input for `tracer summarize`. The output contains summary statistics describing successful TCR reconstruction rates as well information on the cells and which clonal groups they belong to. The output is published to the same S3 bucket.
