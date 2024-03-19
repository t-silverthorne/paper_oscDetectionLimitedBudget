require(dplyr)
require(data.table)
require(simphony)
require(pracma)
featureGroups = data.table(amp = c(0, 1))

B  = 12 


mt = linspace(0,24,B+1)
mt = mt[1:B]
length(mt)


simData = simphony(featureGroups, timepointsType = 'specified',
                   timepoints = mt)
simData %>% {.$abundData} %>% head()

# example for faster shuffling
I=100
data = rnorm(B)
data = replicate(I,{data}) %>% t()
for (ii in c(1:I) ){
  data[ii,]=data[ii,shuffle(data.frame(c(1:B)))]
}

head(data)



x=10
myfun=function(x){
  a=x
  return(a)
}