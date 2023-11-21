library(gurobi)
library(CVXR)
library(dplyr)
N=100
n=10
ncstr=50
nstencil=2
x=Variable(N,boolean=T)

A=matrix(rnorm(N*N),N,N)
A=Constant(A%*%t(A))

L=replicate(ncstr,{sample(c(1,0),N,replace=T)}) %>% t() %>% as.matrix()
Lnum=L
b=rowSums(L)-1
L=Constant(L)
b=Constant(b)

constraints=list(sum(x)==n, sum(pos(L%*%x - b))<=nstencil)

res=solve(Problem(Minimize(quad_form(x,A)),constraints),verbose=T,num_iter=1e6,MIPGapAbs=1e-2,)
x=res$getValue(objet=x)

(Lnum%*%x - rowSums(Lnum) == 0) %>% sum()
