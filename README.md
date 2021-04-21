# Benchmarking DIANE's edge testing procedure

The project [DIANE](https://oceanecsn.github.io/DIANE/) (Dashboard for theInference and Analysis of Networks from Expression data) proposes a novel procedure to sparsify the output of GENIE3 for Gene Regulatory Network inference.

We demonstrate in those scripts that testing regulator-gene pairs by perumtations on random forest feature importance metrics yields to a better precision in the inferred Gene Regulatory Networks than appliying a direct hard-thresholding.

This analysis is carried out on two datasets :

-   A. thaliana's response to heat, salinity and mannitoal stress, validated on the connecTF database.

-   E. coli's response to norfloaxacin in time, validated on the RegulonDB database.

The precision is assessed by the ratio of validated regulatory interactions among predictions over the total predicted interactions with availbale information.

To run the benchmarks, see benchmark_athaliana.R and benchmark_ecoli.R.

To just load the results and plot a final figure, see figure.R.
