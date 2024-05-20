# Overview
There are two components to this repository:

1. Code supporting our manuscript. It can be found in `figs_for_paper/` and `results/` directories.
2. Code that will be packaged as a standalone `R` package. It can be found in the `R/` directory with unit tests in `tests/` and examples in `examples/`.

# Installation of dependencies

## Gurobi
To install Gurobi:

1. Download gurobi
2. Get a license file (free for academic usage)
3. Put the license in /opt/gurobi/ directory (see Windows instructions below)
4. Follow R specific instructions in Gurobi manual to install R package


## Installing all other dependencies:
This package uses `renv` for `R` environment management. To install the remaining dependencies using `renv`, run the following in an `R` console
```R
require(renv)
renv::restore()
```
For more information on using `renv` to manage shared `R` repositories, see the documentation [here](https://rstudio.github.io/renv/articles/collaborating.html).


# Testing package

You can execute the unit tests for this package by running:
```R
require(devtools)
devtools::test()
```