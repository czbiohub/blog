---
layout: post
mathjax: true
title: From transcription factors to reprogramming protocols
date: 2018-10-26
author: Angela Oliveira Pisco
---


In [Tabula Muris](https://www.nature.com/articles/s41586-018-0590-4) we used a random forest model strategy to identify sets of transcription factors using [FACS](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0590-4/MediaObjects/41586_2018_590_MOESM8_ESM.xlsx) and [droplet](https://github.com/czbiohub/tabula-muris/blob/master/23_tf_analysis/rf.model.one.vs.all.transcriptionfactors.droplet.xlsx) data that can potentially be used to inform the design of novel reprogramming protocols.

The [varSelRF](https://cran.r-project.org/web/packages/varSelRF/index.html) package starts with a [random forest](https://en.wikipedia.org/wiki/Random_forest) model to calculate the importance of each gene for defining cell types. Next it uses out-of-bag error as the minimization criterion and carries out variable elimination by successively eliminating the least important variables (with importance as returned from the random forest analysis). The algorithm iteratively fits random forests, at each iteration building a new forest after discarding those variables (genes) with the smallest importance; the selected set of genes is the one that yields the smallest out-of-bag error rate. This leads to the selection of small sets of non-redundant variables.

![Random forest model using transcription factors](/images/reprogramming-direct-diff/rf_tfs_summary.png)

Our confidence in the model comes from comparing the top candidates with the transcription factors currently used to reprogram cell types.

## Contribution of transcription factors to cell identity
 The choice to use transcription factors to pull each cell type apart from the rest was natural, given that using only transcription factors we were able to reconstruct the [full dendrogram of cell identities](https://www-nature-com.ucsf.idm.oclc.org/articles/s41586-018-0590-4/figures/15) with 90% confidence.

![Entanglements](/images/reprogramming-direct-diff/rf_entanglements.png)

Entanglement is a measure of alignment between two dendrograms and the entanglement score ranges from 0 (exact alignment) to 1 (no alignment). To compute the different entanglements we used the dendrogram created from all expressed genes as the reference for comparisons to the dendrograms produced using particular gene ontology cellular functions (transcription factors, cell surface markers, RNA splicing factors).


## Cell surface markers
While the dendrogram obtained when using only cell surface markers is not as good at pulling cell types apart as it is when we use transcription regulators, the cell surface antigens are some of the most commonly used proteins in any cell biology lab. The majority of cell surface markers are molecules within the cell's plasma membrane that are unique to different cell types. We tried the same random forest model to come up with cell sorting panels using the [FACS](https://raw.githubusercontent.com/czbiohub/tabula-muris/blob/master/23_tf_analysis/rf.model.one.vs.all.cellsurfacemarkers.facs.xlsx) and the [droplet](https://raw.githubusercontent.com/czbiohub/tabula-muris/blob/master/23_tf_analysis/rf.model.one.vs.all.cellsurfacemarkers.droplet.xlsx) data in the Tabula Muris.

Despite cell surface markers performance in the above chart, for the majority of the cell types tested we got a small, manageable list of markers that are cell type specific.

## Conclusion
[Induced pluripotent stem cells](https://www.eurostemcell.org/ips-cells-and-reprogramming-turn-any-cell-body-stem-cell) are the holy grail of regenerative medicine; however, the low efficiency of somatic cell reprogramming and the lengthy cell differentiation protocols complicate their potential clinical applications. For this reason, one of the current major goals of single cell research is defining cell identities because there lies the key to the development of better reprogramming protocols. With our machine learning approach we are not only able to provide candidate transcription factors for novel reprogramming protocols, but we can already take a step further and suggest cell type validation sets, as our model can also predict cell surface markers that be used either for FACS sorting or antibody staining.

While this is still very experimental, the results match previously known findings, as for example the combination of *Cd36* and *Cav1* being enough to sort fat endothelial cells. We are interested in knowing whether our model predictions are valid in cell types not previously reported, so if you take this further and experimental validate the results by either testing the transcription factors or the cell surface markers candidates reach out!

*[@aopisco](https://github.com/aopisco) on the behalf of the Tabula Muris consortium*
