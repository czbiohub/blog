---
layout: post
mathjax: true
codehide: true
title: How can scRNA-seq inform the development of new reprogramming protocols and sorting markers pannels?
date: 2018-10-26
author: Angela Oliveira Pisco
---



## From transcription factors to reprogramming protocols
[Induced pluripotent stem cells](https://www.eurostemcell.org/ips-cells-and-reprogramming-turn-any-cell-body-stem-cell) are the holy grail of regenerative medicine; however, the low efficiency of somatic cell reprogramming and the lengthy cell differentiation protocols are hardening their potential clinical applications. For this reason, one of the current major goals of single cell research is defining cell identities because there lies the key to the development of better reprogramming protocols. Critical to such work is the understanding of the underlying regulatory networks. While a handful of molecular functions are usually associated with cell regulation, the most commonly used are transcription factors or transcription regulators in general. In [Tabula muris](https://www.nature.com/articles/s41586-018-0590-4) we investigated how transcription factors are able to reconstruct known cell ontology relationships between cell populations. Using a random forest model strategy we identified sets of transcription factors using [facs](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0590-4/MediaObjects/41586_2018_590_MOESM8_ESM.xlsx) and [droplet](/files/rf.model.one.vs.all.cellsurfacemarkers.droplet.xlsx) data that can potentially be used to inform the design of novel reprogramming protocols.


![Random forest model using transcription factors](/images/reprogramming-direct-diff/rf_tfs_summary.png)


## Contribution of transcription factors to cell identity
Our confidence in the model comes from empirical validation, when comparing the top 10 candidates with the transcription factors currently used to reprogram some of the cell types. The choice of using transcription factors to pull each cell type a part from the rest was natural given that using only transcription factors we were able to reconstruct the full dendrogram of cell identities with 90\% confidence

![Tanglegram transcription factors](/images/reprogramming-direct-diff/rf_tfs_entanglements.png)

in contrast with many other chosen functional groups.

![Entanglements](/images/reprogramming-direct-diff/rf_entanglements.png)

## Cell surface markers
While it is of note the inability of cell surface markers to pull cell types apart, the cell surface antigens are among the most ubiquitous class of proteins used on a daily basis in any cell biology lab. The majority of cell surface markers are molecules within cell's plasma membrane that are unique to different cell types. We sought to try using the random forest model used for identifying transcription factors associated with individual cell identities to come up with cell sorting panels using the [facs](/files/rf.model.one.vs.all.cellsurfacemarkers.facs.xlsx) and the [droplet](/files/rf.model.one.vs.all.cellsurfacemarkers.droplet.xlsx) data in Tabula muris.


## Seeking collaborations
While this model is still very experimental, the results match previously known findings, as for example the combination of *Cd36* and *Cav1* being enough to sort fat endothelial cells. We would like to take this further and experimental validate the results so if you are willing to test the transcription factors or cell surface markers candidates reach out!

*[@aopisco](https://github.com/aopisco) on the behalf of the Tabula muris consortium*
