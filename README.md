# Benchmarking DIANE's edge testing procedure


The project **DIANE** (Dashboard for theInference and Analysis of Networks from Expression data) proposes a novel procedure to sparsify the output of the GENIE3 for Gene Regulatory Network inference.

We demonstrate in those scripts that testing regulator-gene pairs by perumtations on random forest feature importance metrics yields to a better precision in the inferred Gene Regulatory Networks.

This is done on two datasets : A. thaliana's response to heat, and E. coli's response to norfloxacin. The precision is assessed on available known regulatory interactions in the databases connecTF and RegulonDB respectively.

To run the benchmarks, see benchmark_athaliana.R and benchmark_ecoli.R.
To just load the results and plot a final figure, see figure.R.