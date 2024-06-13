# scRNA-seq Integrated Cancer Attractor Analysis using Python

[![Coverage Status](https://coveralls.io/repos/github/marcosgvjunior/scRNA-seq-Integrated-Cancer-Attractor-Analysis-using-Python/badge.svg)](https://coveralls.io/github/marcosgvjunior/scRNA-seq-Integrated-Cancer-Attractor-Analysis-using-Python)

## Overview
Repository to progressively translate the Wolfram Mathematica code from the original repository ([`Biomarker-Guided scRNA-Seq Cancer Attractor Analysis`](https://github.com/marcosgvjunior/Biomarker-Guided-scRNA-Seq-Cancer-Attractor-Analysis)) to Python. It also aims to include analyses from previous studies ([article](https://doi.org/10.3390/ijms25094894)) and additional features.

## Features Already Implemented
- Construct and visualize gene regulatory networks (GRNs).
- Generate marker gene combinations and visualize dispersion.
- Create histograms and t-SNE plots for dimensionality reduction.
- Apply clustering algorithms (KMeans, DBSCAN) on reduced/markers dimensions and analyze clusters.

## Features to Implement
- Investigate biomarkers' potential for finding cancer attractors.
- Perform stochastic simulations to refine attractor identification.
- Identify potential clusters' multistability as sources of cancer recurrence.
- Others: Additional analyses and features from previous and upcoming studies.

## Input Files

1. Gene expression data: 
   - "inputs/datapoints_seurat_BT_ALL.xlsx"
   
   This dataset contains gene expression data. The data can be obtained from the [link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84465).

2. Activation and inhibition adjacency matrices: 
   - "inputs/act_BT_ALL_seurat.xlsx"
   - "inputs/sup_BT_ALL_seurat.xlsx"
   
   These are the activation and inhibition adjacency matrices corresponding to the GRN dynamic model.

## Usage

1. **Data**: Input scRNA-seq data in a normalized count matrix format with adjacency matrices corresponding to the GRNs’ activation and inhibition interactions.
2. **Dimensionality reduction and clustering**: Select the desired biomarker dimensions to proceed with the analysis.

## Dependencies
- Python 3.x
- Jupyter Notebook
- NumPy
- Pandas
- Matplotlib
- SciPy
- Scikit-learn
- NetworkX
- SymPy
- TQDM
- Joblib

## Requirements

All necessary packages are pre-installed in Google Colab.

To view and run the notebook:

1. Clone the repository to your local machine.
2. Upload the repository to Google Colab.
3. Run the cells as needed.

## License

This project is licensed under the [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License (CC BY-NC-SA 4.0)](LICENSE.txt).

### License Summary

This license allows others to remix, adapt, and build upon this work non-commercially, as long as they credit the original author and license their new creations under the same terms. For more details, please visit the [license page](https://creativecommons.org/licenses/by-nc-sa/4.0/).

## Contribution
Any contributions are welcome.

## Authors
- [Marcos Vieira](https://github.com/marcosgvjunior)
