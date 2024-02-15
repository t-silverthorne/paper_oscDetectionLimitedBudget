# count partitions using recurrence relation
require(devtools)
devtools::load_all('.')

# partition counting
sa_enumpar(5) # expect 1 1 2 3 5 7
sa_enumpar(50) # larger example

# generate sampling probabilities for partition
sa_eulermat(30)

