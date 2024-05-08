# SpcvxBic
Identifies a bicluster of a data matrix, clustering observations and features simultaneously, and apply sparsity to features. 'SpcvxBic' can be used for high-dimensional data, identifies uninformative features. The method is based on Bi-ADMM algorithm. For a detailed introduction to this algorithm, please refer to the author's related papers.

# Note
Download the compressed file with the extension ".gz" to locally install the SpcvxBic package in RStudio.
main function is "SCB_ADMM_speed.R". Please use the corresponding Python file and place it in the working directory to smoothly call Python functions and accelerate program computations. "SCB_ADMM_speed_WS.R" can use previous grid search result {A,v,g,z} to accelerate calculate speed.

# Example
Full code of simulation were provided with supplyment file of paper "sparse convex bicluster".
