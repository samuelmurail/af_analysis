---
title: 'Af-analysis: a Python package for Alphafold analysis'
tags:
  - Python
  - Alphafold
  - Protein Structure
  - Structural bioinformatics
authors:
  - name: Alaa Reguei
    affiliation: 1
  - name: Samuel Murail
    orcid: "0000-0002-3842-5333"
    corresponding: true
    affiliation: "1, 2"

affiliations:
 - name: Université Paris Cité, Inserm, CNRS, BFA, F-75013 Paris, France
   index: 1
 - name: Ressource Parisienne en Bioinformatique Structurale (RPBS), F-75013 Paris, France
   index: 2

date: 8 January 2025
bibliography: paper.bib
---

# Summary

The publication of AlphaFold 2 [@jumper_highly_2021] has significantly advanced the field of protein structure prediction. The prediction of protein structures has long been a central challenge in the field of structural bioinformatics, with the ultimate goal of elucidating the relationship between protein structure and function. Accurate prediction of protein structure is essential for a number of applications, including drug discovery, protein engineering, and the study of protein-protein interactions. AlphaFold, which employs a deep learning-based approach, has demonstrated unprecedented accuracy in protein structure prediction, outperforming other contemporary methods. In this paper, we present `af-analysis`, a Python package that provides tools for the analysis of AlphaFold results. The `af-analysis` library has been designed to facilitate the analysis of many different protein structures predicted by AlphaFold and its derivatives. It provides functions for comparing predicted structures with experimental structures, visualising predicted structures, and calculating structural quality metrics.

# Statement of need

With the release of AlphaFold 2 [@jumper_highly_2021] in 2021, the scientific community has achieved an unprecedented level of accuracy in predicting protein structures. Derivatives of AlphaFold 2, namely ColabFold [@mirdita2022colabfold], AlphaFold Multimer [@Evans2021.10.04.463034], AlphaFold 3 [@abramson_accurate_2024] and its re-implementations such as Boltz-1 [@wohlwend_boltz-1_2024] and Chai-1 [@discovery_chai-1_2024] have been developed to predict the structure of protein complexes, setting a new standard for protein-protein and protein-peptide docking.

Analysis of AlphaFold results is a crucial step in the process of utilising these predictions for scientific research. The AlphaFold software provides several excellent quality metrics that offer valuable information about the accuracy of the predicted structures. Among these scores, the predicted local distance difference test (pLDDT) is a per-residue measure of local confidence, as the predicted aligned error (PAE) provides confidence over the relative position of two residues within the predicted structure. To analyze these results, the AlphaBridge webserver [@alvarez-salmoral_alphabridge_2024] and the PICKLUSTER plugin [@genz_pickluster_2023] for the UCSF ChimeraX visualisation software were developed to characterize the different interfaces within protein complexes, and extract their respective scores.

Although these tools are very practical, Bjorn Wallner has shown that calculating 5 or 25 basic AlphaFold models may not be enough, it is sometimes necessary to generate thousands of models to obtain a few high quality models, leading to the AlphaFold derivative, AFsample [@10.1093/bioinformatics/btad573]. Massive sampling altogether with multiple software usage (AFsample and ColabFold), weights and parameters has been integrated into the MassiveFold software [@raouraoua_massivefold_2024] and has shown performance approaching the accuracy of AlphaFold 3.


The subsequent analysis of hundreds to thousands of models can prove to be a tedious and meticulous process, as dealing with thousands of models and different output formats can be time consuming.
Furthermore, if the quality metrics produced by AlphaFold are good, additional metrics have been developed to assess the quality of the models. These include pdockq [@bryant2022improved], pdockq2 [@10.1093/bioinformatics/btad424], and LIS score [@Kim2024.02.19.580970]. All of the these metrics have to be calculated from different scripts. Another point to consider is the diversity of the models. As shown in AFsample, it is sometimes necessary to compute up to tens of thousands of models and then cluster them in order to select the best ones. The `af-analysis` library has been developed to facilitate the analysis of sets of model structures and associated metrics. The library is based on the `pandas` library and is able to import AlphaFold 2 and 3, ColabFold, Boltz-1 and Chai-1 prediction directory as a `pandas` DataFrame. The library provides a number of functions to add further metrics to the DataFrame, compare models with experimental structures, visualise models, cluster models and select the best models.

# References