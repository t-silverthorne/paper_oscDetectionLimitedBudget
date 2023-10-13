# testing regularization
library(stats)
library(dplyr)
library(ggplot2)
library(gurobi)
library(CVXR)
library(dplyr)
library(expm)
library(patchwork)
library(pracma)
library(data.table)
source('utils/powerUtils.R')

# define problem
N=20
param=list(Amp         = 1.5,
           freq        = 4.4,
           acro        = 1,
           method      = 'min',
           regStrength = 1)
t=(0:N)/N
t=t[1:(length(t)-1)]
t0=t

fmin = 1 
fmax = 24
df   = .1 
freqlist = seq(from=fmin,to=fmax,by=df) %>% as.list()

vis_colors <- c("poly" = "blue", "seq" = "red", 'unif'='black')



# iterate over partitions to find optimal solutions
N1list = seq(from=4,to=floor(N/2)) %>% as.list()
varyN1df= lapply(N1list,function(N1){
  N2=N-N1
  tvec1 = (0:N1)/N1 
  tvec2 = (0:N2)/N2 
  tvec1 = tvec1[1:(length(tvec1)-1)]
  tvec2 = tvec2[1:(length(tvec2)-1)]
  
  res_delta = pracma::fminbnd(function(delta){-Jfun_delta(delta,tvec1,tvec2,freqlist,param)},
                    a=0,
                    b=max(diff(tvec2)))
  res_tau   = pracma::fminbnd(function(tau){-Jfun_tau(tau,tvec1,tvec2,freqlist,param)},
                    a=0,
                    b=1)

   return(list(delta_opt  = res_delta$xmin,
               delta_fmin = res_delta$fmin,
               tau_opt    = res_tau$xmin,
               tau_fmin   = res_tau$fmin,
               N1         = N1))
}) %>% rbindlist()

# extract optimal solutions
delta_opt = varyN1df[which.min(varyN1df$delta_fmin),]$delta_opt
N1poly   = varyN1df[which.min(varyN1df$delta_fmin),]$N1
N2poly   = N - N1poly
tvec1_poly = (0:N1poly)/N1poly
tvec2_poly = (0:N2poly)/N2poly
tvec1_poly = tvec1_poly[1:(length(tvec1_poly)-1)]
tvec2_poly = tvec2_poly[1:(length(tvec2_poly)-1)]

tau_opt   = varyN1df[which.min(varyN1df$tau_fmin),]$tau_opt
N1seq     = varyN1df[which.min(varyN1df$tau_fmin),]$N1
N2seq     = N - N1seq
tvec1_seq = (0:N1seq)/N1seq
tvec2_seq = (0:N2seq)/N2seq
tvec1_seq = tvec1_seq[1:(length(tvec1_seq)-1)]
tvec2_seq = tvec2_seq[1:(length(tvec2_seq)-1)]

tvec_poly = construct_poly_design(delta_opt,tvec1_poly,tvec2_poly)
tvec_seq  = construct_sequential_design(tau_opt,tvec1_seq,tvec2_seq)


# extract optimal solutions
delta_opt = varyN1df[which.min(varyN1df$delta_fmin),]$delta_opt
N1poly   = varyN1df[which.min(varyN1df$delta_fmin),]$N1
N2poly   = N - N1poly
tvec1_poly = (0:N1poly)/N1poly
tvec2_poly = (0:N2poly)/N2poly
tvec1_poly = tvec1_poly[1:(length(tvec1_poly)-1)]
tvec2_poly = tvec2_poly[1:(length(tvec2_poly)-1)]

tau_opt   = varyN1df[which.min(varyN1df$tau_fmin),]$tau_opt
N1seq     = varyN1df[which.min(varyN1df$tau_fmin),]$N1
N2seq     = N - N1seq
tvec1_seq = (0:N1seq)/N1seq
tvec2_seq = (0:N2seq)/N2seq
tvec1_seq = tvec1_seq[1:(length(tvec1_seq)-1)]
tvec2_seq = tvec2_seq[1:(length(tvec2_seq)-1)]

tvec_poly = construct_poly_design(delta_opt,tvec1_poly,tvec2_poly)
tvec_seq  = construct_sequential_design(tau_opt,tvec1_seq,tvec2_seq)

# Generate first panel: a plt of the designs
tdf = list(data.frame(time=tvec_poly,value=2,type='poly'),
     data.frame(time=(tvec1_poly+delta_opt)%%1,value=1,type='unif'),
     data.frame(time=tvec2_poly,value=3,type='unif')) %>% rbindlist()
polyA = tdf %>% ggplot(aes(x=time,y=value,color=type))+geom_point()+scale_color_manual(values=vis_colors)

tau=res_tau$xmin  
tdf = list(data.frame(time=tvec_seq,value=2,type='seq'),
     data.frame(time=tvec1_seq*tau_opt,value=1,type='unif'),
     data.frame(time=tau_opt+tvec2_seq*(1-tau_opt),value=3,type='unif')) %>% rbindlist()
seqA  = tdf %>% ggplot(aes(x=time,y=value,color=type))+geom_point()+scale_color_manual(values=vis_colors)

tdf = list(data.frame(time=t0,value=2,type='unif')) %>% rbindlist()
unifA  = tdf %>% ggplot(aes(x=time,y=value,color=type))+geom_point()+scale_color_manual(values=vis_colors)

polyA  = polyA+theme(legend.position = "none")
seqA   = seqA+theme(legend.position = "none")
unifA  = unifA+theme(legend.position = "none")
panelA = ((polyA+ggtitle('Design'))/seqA/unifA)

############# move this eventually
df_plotC = 1e-2

# Generate third panel: plot of eigenvalue
freqs=seq(from=fmin,to=fmax,by=df_plotC)

hmatC = freqs %>%as.list() %>%  lapply(function(frq){
    param$freq = frq
    eigUnif = getMinEig(t0,param)
    eigSeq  = getMinEig(tvec_seq,param)
    eigPoly = getMinEig(tvec_poly,param)
    df1 = data.frame(freq=frq,minEig=eigUnif,sched='unif')
    df2 = data.frame(freq=frq,minEig=eigSeq,sched='seq')
    df3 = data.frame(freq=frq,minEig=eigPoly,sched='poly')
    return(rbind(df1,df2,df3))
}) %>% rbindlist()

panelC <- hmatC %>% ggplot(aes(x=freq,y=minEig,color=sched,group=sched))+geom_line()+scale_color_manual(values=vis_colors)+facet_wrap(~sched,nrow=3)


panelA |  panelC+ggtitle('Minimal eigenvalue')
