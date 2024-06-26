[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
# Fermat Distances: Metric Approximation, Spectral Convergence, and Clustering Algorithms

This repository contains all the code necessary to replicate the experiments contained in the paper *Fermat Distances: Metric Approximation, Spectral Convergence, and Clustering Algorithms* ([arxiv link](https://arxiv.org/pdf/2307.05750)) by [Nicolás García Trillos](https://www.nicolasgarciat.com/), [Anna Little](https://www.anna-little.com/), [Daniel McKenzie](https://danielmckenzie.github.io/), and [James M. Murphy](https://jmurphy.math.tufts.edu/). 

## Installation
You will need to compile the C++ files in ```Auxiliary/PowerWeightedShortestPaths/Modified_Dijkstra/Priority_Queue``` before using this code. For this you will need a compiler, which is platform-specific. See [here](https://www.mathworks.com/support/requirements/supported-compilers.html) for further details.

## Quickstart
For a quick demonstration, run either file in ```ImageSegmentationExperiments```.

## Experiments
 - Image segmentation experiments (Figures 1, 4) can be generated by running ```/ImageSegmentationExperiments/BlueSky_Script.m``` or ```/ImageSegmentationExperiments/Moon_Script.m```.
 - Fermat geodesic images (Figure 2) can be generated by running ```/Fermat_Geodesics/Plot_Geodesic.m```
 - Clustering experiments on elongated data with a density gap (Figure 3) can be generated via ```/GeometricVersusDensityClustering/GeometricVersusDensity_Script.m```
 - Eigenvalue comparison experiments (Figures 5-12) are generated via ```/ComparisonExperiments/BulkPlotsScript.m```, choosing which of the two data sets to use.

## Auxiliary code
The subdirectory ```/Auxiliary/``` contains a minimal version of the ```PowerWeightedShortertPaths``` package developed by Daniel Mckenzie and [Steven Damelin](https://scholar.google.com/citations?hl=en&user=nVqG2rwAAAAJ&view_op=list_works&sortby=pubdate).See *Power weighted shortest paths for clustering Euclidean data.* ([arxiv link](https://arxiv.org/pdf/1905.13345)) for more information.  It also contains an open source implementation of the Hungarian algorithm, courtesy of Yi Cao.  

## Citation
If you find this code useful, please cite our paper:
```
@article{trillos2023fermat,
  title={Fermat Distances: Metric Approximation, Spectral Convergence, and Clustering Algorithms},
  author={Trillos, Nicol{\'a}s Garc{\'\i}a and Little, Anna and McKenzie, Daniel and Murphy, James M},
  journal={arXiv preprint arXiv:2307.05750},
  year={2023}
}
```

