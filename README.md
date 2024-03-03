# Overview
There are two components to this repository:

1. Supporting code for our manuscript can be found in `figs_for_paper/` and `results/` directories.
2. Code that will eventually be packaged as a standalone `R` package, found in the `R/` directory with unit tests in `tests/`.

Closer to the release of the paper, items (1) and (2) will be split into separate repositories. 

# Installation of non-standard packages

Two non-standard R packages are used in this repository and their installation instructions are given below. 

Before installation, the user should download Rstudio and open the `oscDetectPaper.Rproj` file located in the root directory of this repo. This will ensure that all dependencies are installed locally and will not affect any pre-existing R environments.

## CVXR
We require a specific branch of the CVXR repo. It can be installed by opening an R console running the following: 
```R
# Make sure you have devtools
require(devtools)

# Specify the GitHub repository
github_repo <- "cvxgrp/CVXR"

# Specify the branch name 
github_branch <- "WIP_gurobi_solver_status"

# Install the package from the specific branch
devtools::install_github(repo = github_repo, ref = github_branch)
```

To test your installation, run:
```R
require(CVXR)
```

## Gurobi

Solving mixed-integer programming problems in CVXR requires access to the commercial Gurobi solver. To install Gurobi and link it to CVXR:

1. Download gurobi
2. Get license file
3. Put license in /opt/gurobi/ directory (see Windows instructions below)
4. Follow R specific instructions in Gurobi manual to install R package
5. Follow OS specific instructions to setup R dynamic library loading 

### Mac and Linux
On Linux, the .profile file should be updated to read
```bash
export GUROBI_HOME="/opt/gurobi1003/linux64"
export PATH="${PATH}:${GUROBI_HOME}/bin"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"
```
the first line should be modified appropriately for users on a Mac operating system to the corresponding path where Gurobi is installed.

### Windows
No modification of your path is necessary provided that your license file is placed in the following location. 
```bash
C:/gurobi/gurobi.lic
```

## Testing dynamic library loading

To make sure the previous step worked (i.e. Gurobi can be accessed as a backend solver for CVXR), open an R console and run the following.

```R
source('examples/ex_cvxr.Rmd')
```

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