library(gurobi)
library(CVXR)
library(dplyr)
library(pracma)
N=50
n=10
ncstr=80
nstencil=2
x=Variable(N,boolean=T)

A=matrix(rnorm(N*N),N,N)
A=Constant(A%*%t(A))

L=replicate(ncstr,{sample(c(1,0),N,replace=T)}) %>% t() %>% as.matrix()
Lnum=L
b=rowSums(L)-1
bnum=b
L=Constant(L)
b=Constant(b)

c1 = sample(c(1,0),N,replace=T)
c2 = sample(c(1,0),N,replace=T)
max_elemwise(max_entries(conv(Constant(c1),x)),max_entries(conv(Constant(c2),x)))

constraints=list(sum(x)==n, 
                 sum(pos(L%*%x - b))<=Constant(nstencil),
                 sum(neg(Lnum%*%x-rowSums(Lnum)))<=5e3)

#constraints=list(sum(x)==n, 
#                 sum(pos(L%*%x - b))<=Constant(nstencil),
#                 sum(neg(Lnum%*%(x-t(Lnum)%*%pos(Lnum%*%x - b))-rowSums(Lnum)))<=5e3)

res=solve(Problem(Minimize(quad_form(x,A)),constraints),verbose=T,num_iter=1e6,MIPGapAbs=1e-2,)
x=res$getValue(objet=x)

