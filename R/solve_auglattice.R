solve_auglattice=function(dvar0,freqs,control,...){
  N1_init     = dvar0[['N1']]
  N2_init     = dvar0[['N2']]
  shift2_init = dvar0[['shift2']]
  scale2_init = dvar0[['scale2']]
 
  x0 = c(N1_init,N2_init,shift2_init,scale2_init)
  cfun = function(x){-costfun_auglattice(N1     = x[1],
                                         N2     = x[2],
                                         shift2 = x[3],
                                         scale2 = x[4],
                                         freqs  = freqs,
                                         ...)}
  tfun = function(x){tfun_auglattice(N1     = x[1],
                                     N2     = x[2],
                                     shift2 = x[3],
                                     scale2 = x[4],
                                     freqs  = freqs,
                                     control=control)}
   xopt=stats::optim(x0,fn=cfun,gr=tfun,
                  method='SANN',
                  control=list(trace  = control$trace,
                               REPORT = control$REPORT,
                               maxit  = control$maxit))
   return(xopt)
}