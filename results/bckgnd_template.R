library(stringr)
library(lubridate)
library(job)
library(gurobi)
library(CVXR)
gc()
library(devtools)
load_all('.')
test=F

tstamp <- now() %>% toString() %>% str_replace(' ','___')

##### no lattice constraints

# simulated annealing runs
if (test){
  ens_size = 10 
} else {
  ens_size = 30 
}

tstamp = now() %>% toString() %>% str_replace(' ','___')
outdir =paste0('results/output/batch_',tstamp,'/')
dir.create(outdir)

for (Nmeas_now in c(20:20)){
  for (solver in c('sa')){
    if (solver=='sa'){
      load_all('.')
      opts=make_default_opts(prob_size='medium',solver_type='simulanneal')
      opts$lattice_cstr = 'sa_lattice'
      if (opts$lattice_cstr=='sa_lattice'){
        opts$min_active_lats = 2
        opts$max_active_lats = 2
        opts$min_lat         = 4
        opts$max_lat         = 10
      }
      
      opts$Nmeas = Nmeas_now
      opts$Nfine        = 12*24 # 5 minute intervals 
      opts$fmin         = 1
      opts$fmax         = 24
      opts$costfun_type = 'Linfty'
      opts$verbose      = F
      if (test){
        opts$Nfreq      = 2^3
        opts$num_iter   = 10
      }else{
        opts$Nfreq      = 2^10
        opts$num_iter   = 1e2
      }
      
            
      if (opts$lattice_cstr=='sa_lattice'){
        sa_sols=as.list(replicate(ens_size,{NULL})) 
      }else if(opts$lattice_cstr=='none'){
        sa_sols = replicate(ens_size*opts$Nmeas, {NaN})
        sa_sols = matrix(sa_sols,nrow=ens_size,ncol=opts$Nmeas)
      }else{
        stop('unrec')
      }
      
      
      Aquad=make_quadmats(opts)
      tau          = c(0:opts$Nfine)/opts$Nfine       
      tau          = tau[1:(length(tau)-1)] 
      
      for (ii in c(1:ens_size)){
        print(toString(ii))
        xopt         = run_sa_power(Aquad,opts)
      
        if (opts$lattice_cstr=='none'){
          mt_opt_sa    = tau[xopt$xval==1]
          sa_sols[ii,] = mt_opt_sa
        }else if(opts$lattice_cstr=='sa_lattice'){
          sa_sols[[ii]] = xopt
        }else{
          stop('unrecognized')
        }
      }
      tstamp = now() %>% toString() %>% str_replace(' ','___')
      fname  = paste0(outdir,tstamp,'_simulanneal_Nmeas_',toString(Nmeas_now))
      saveRDS(sa_sols,paste0(fname,'_solution.RDS'))
      saveRDS(opts,paste0(fname,'_opts.RDS'))
      
    }else if(solver=='cvxr'){
      opts=make_default_opts(prob_size='medium',solver_type='cvxr')
      opts$Nmeas = Nmeas_now
      opts$Nfine        = 12*24 # 5 minute intervals 
      opts$fmin         = 1
      opts$fmax         = 24
      opts$costfun_type = 'Linfty'
      opts$verbose      = T
      if (test){
        opts$Nfreq      = 2^3
        opts$time_limit = 5
      }else{
        opts$Nfreq      = 2^8
        opts$time_limit = 60*30 # 1/2 hr
      }
      opts$MIPFocus=1
      opts$Presolve=2
      x     = make_variable(opts)
      csts  = make_constraints(x,NULL,NULL,opts)
      
      Aquad = make_quadmats(opts)
      prob  = make_problem(x,Aquad,csts,opts)
      start=Sys.time()
      xopt = CVXR::solve(prob,verbose=opts$verbose,num_iter=1e9,
                         TimeLimit=opts$time_limit,
                         MIPGapAbs=opts$MIPGapAbs,
                         MIPFocus=opts$MIPFocus,
                         Presolve=opts$Presolve
                          )
      tstamp = now() %>% toString() %>% str_replace(' ','___')
      fname  = paste0(outdir,tstamp,'_cvxr_Nmeas_',toString(Nmeas_now))
      # saveRDS(xopt[[1]],paste0(fname,'_solution.RDS'))
      saveRDS(opts,paste0(fname,'_opts.RDS'))
      
    }
  }
}
note='completed'
write(note,paste0(outdir,'completed.txt'))

# cvxr runs

