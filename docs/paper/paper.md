---
title: 'Af-analysis: a Python package for Alphafold analysis'
tags:
  - Python
  - Alphafold
  - Protein Structure
  - Structural bioinformatics
authors:
  - name: Alaa Reguei
    affiliation: "1"
  - name: Samuel Murail
    orcid: "0000-0002-3842-5333"
    affiliation: "1"

affiliations:
 - name: "Université Paris Cité, Inserm U1133, CNRS UMR 8251, Paris, France"
   index: 1

date: 11 July 2024
bibliography: paper.bib
---

# Summary

The publication of AlphaFold 2[@jumper_highly_2021] has significantly advanced the field of protein structure prediction. The Prediction of protein structures has long constitued a central challenge in the field of structural bioinformatics, with the ultimate objective of elucidating the relationship between protein structure and function. The accurate prediction of protein structures is of great importance for a number of  applications, including drug discovery, protein engineering, and the investigation of protein-protein interactions. AlphaFold 2, which employs a deep learning-based approach, has demonstrated unprecedented accuracy in protein structure prediction, outperforming other contemporary methods. In this paper, we present `af-analysis`, a Python package that provides tools for the analysis of Alphafold 2 results. `af-analysis` has been designed to facilitate the analysis of protein structures predicted by Alphafold. It provides functions for the comparison of predicted structures with experimental structures, the visualisation of predicted structures, and the calculation of structural quality metrics.

# Statement of need

The publication of AlphaFold 2[@jumper_highly_2021] in 2021, has enabled the scientific community to achieve a previously unattainable level of accuracy in predicting protein structures. Derivatives of AlphaFold 2, namely Colabfold[@mirdita2022colabfold] and AlphaFold Multimer[@Evans2021.10.04.463034] have been developed to predict the structure of protein complexes, thereby establishing a new standard for protein-protein and protein-peptide docking. However as demonstrated by Björn Wallner, the generation of thousands of models is sometimes necessary to obtain a few high-quality models[@10.1093/bioinformatics/btad573]. Subsequent analysis of the results may prove to be a laborious and meticulous process. Furthermore if quality metrics produced by AlphaFold are good, additional metrics have been developed to assess the quality of the models. This include pdockq[@bryant2022improved], pdockq2[@10.1093/bioinformatics/btad424], and LIS score[@Kim2024.02.19.580970]. All of the aforementioned metrics must be calculated from different scripts. A further point to consider is the diversity of the models. As illustrated in AFsample[@10.1093/bioinformatics/btad573], there are occasions when it is necessary to compute up to ten thousand models, and then cluster them in order to select the most appropriate models. The `af-analysis` has been developed with the objective of facilitating the analysis of sets of model structures and associated metrics. The library has been constructed on the `pandas` library and is capable of importing an alphafold or colabfold prediction directory as a `pandas` DataFrame. The library provides a range of functions that enable addition of further metrics to the DataFrame, the comparison of models with experimental structures, the visualisation of the models, the clustering of models, and the selection of the best models.

# References