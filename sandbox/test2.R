library(devtools)
load_all()
lat1=c(1:10)/10-1/10
lat2=c(1:5)/5 -1/5
control=NULL
control$tfun_choice  = 'default'
control$mean_scale1  = 0 
control$mean_scale2  = 0 
control$mean_shift2  = 0 
control$sd_scale1    = .2
control$sd_scale2    = .1
control$sd_shift2    = .1

shift1=0
shift2=.5
scale1=.5
scale2=.5

tfun_2lattice_svdpower(shift1=shift1,shift2=shift2,scale1=scale1,scale2=scale2,lat1,lat2,control)
