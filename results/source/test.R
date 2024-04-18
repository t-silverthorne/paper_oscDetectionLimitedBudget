require(dplyr)
require(gurobi)
source('results/source/powerChord.R')
powerChord(38,drts=Inf,w_reg=1e-2,tlim=20)