# Bayesian approaches to designing replication studies

This repository contains code and data related to the paper

Pawel, S., Consonni, G., Held, L. (2023). Bayesian approaches to designing
replication studies. Psychological Methods.
DOI:[10.1037/met0000604](http://doi.org/10.1037/met0000604). Preprint: <https://arxiv.org/abs/2211.02552>

## Reproducing the results

We offer two ways to reproduce the results

### 1. Reproduction with local computational environment (requires R and LaTeX)

First install the required R packages by running in a shell from the root
directory of the repository

``` sh
## packages from CRAN
R -e 'install.packages(read.delim("CRANpackages.txt", header = FALSE)[,1])'
## requires remotes package, also available on CRAN
R -e 'remotes::install_gitlab(repo = "samuel.pawel/BayesRep", subdir = "pkg",
                              host = "gitlab.uzh.ch")'
```

Then run

``` sh
cd paper
make pdf
```

this should reproduce all analyses and output the file `batdrs.pdf` in the paper
directory.

Although our analysis depends on only few dependencies, this approach may lead
to different results (or not even compile successfully) in the future if R or an
R package dependency changes. The R and R package versions which were used in
our analysis can be seen in the output of the sessionInfo command at the bottom
of the manuscript in the snapshot of the GitHub repository at the time of
submission.

### 2. Reproduction within Docker container (requires Docker with root rights)

Run in a shell from the root directory of the repository

``` sh
make drunpdf
```

this should output the file `batdrs.pdf` in the paper directory. The Docker
approach takes a bit longer but reruns our analyses in a Docker container which
encapsulates the computational environment (R and R package versions) that was
used in the original analysis. 
