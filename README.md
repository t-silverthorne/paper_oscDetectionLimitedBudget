# Overview
This private repository contains code for optimizing rhythm detection experiments with frequency uncertainty. There are two main components to this repository:

1. Code for supporting our manuscript, found mainly in figs_for_paper/ and results/ directories.
2. Code that will eventually be packaged as a standalone R package, found in the R/ directory with unit tests in tests/.

Closer to the release of the paper, items (1) and (2) will likely be split into separate repositories. 


# Installation of non-standard packages
Two non-standard R packages are used in this repository and their installation instructions are given below.

## CVXR
We require a specific branch of the CVXR repo. It can be installed by starting an R session and running the following: 
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
3. Put license in /opt/gurobi/ directory
4. Follow R specific instructions in Gurobi manual to install R package
5. Update .profile so that dynamic library loading works with Gurobi. 

On linux, the .profile file will read
```bash
export GUROBI_HOME="/opt/gurobi1003/linux64"
export PATH="${PATH}:${GUROBI_HOME}/bin"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"
```

the first line should be modified appropriately on Mac or Windows operating systems.


To make sure that Gurobi can be accessed as a backend solver for CVXR, open an R terminal and run the following.
```R
source('examples/ex_cvxr.Rmd')
```

## Testing all other dependencies:
In Rstudio, you can execute the unit tests for this package by running

  require(devtools)
  devtools::test()

If all dependencies are working, the tests should all pass.
