clear all
clf 
x=rand*10;
lambda=2
Nmeas=randsample(10:100,1)
dlv=.001
lv=0:dlv:40
plot(lv(1:end-1),abs(diff(ncfcdf(x,2,Nmeas-3,lv))/dlv))

d1=2;
d2=Nmeas-3;

Iterm=@(x,d1,d2,k) betainc(d1*x/(d2+d1*x),d1/2+k,d2/2)
hold on
%yline(abs(0.5*Iterm(x,d1,d2,0)))
yline(1/Nmeas^0.5)
ylim([0,1])
%%
syms d1 d2 x t
assume(x>0)
assume(d1==2)
assume(d2>0)
betaF=gamma(d1/2)*gamma(d2/2)/gamma(d1/2+d2/2);
limit(simplify(int((1-t)^(d2/2-1),0,d1*x/(d2+d1*x)))/betaF,d2,Inf)

%%
n=10

Nfreq = 2^8;
Nacro = 2^8;
fmin  = 1;
fmax  = 24;

freq  = linspace(fmin,fmax,Nfreq);

dfreq=(fmax-fmin)/(Nfreq-1);
dacro = 2*pi/Nacro;
%%
acro  = linspace(0,2*pi,Nacro+1);
acro  = acro(1:end-1);

[F,A] = meshgrid(freq,acro);
Amp=5

Nvals = 10.^linspace(1,5,5)
eps=1e-3
bvec=[]
for N=Nvals
    %nexttile
    mt=rand(1,N);
    rand_pert = eps*rand(1,N);
    bvec(end+1)=norm(get_beta(mt+rand_pert,F,A,Amp,dfreq,dacro)-get_beta(mt,F,A,Amp,dfreq,dacro))/norm(rand_pert);
end
bvec

function beta_val=get_beta(mt,F,A,Amp,dfreq,dacro)
    alpha=.05;
    Nmeas = length(mt);
    mt       = reshape(mt,1,1,length(mt));
    csq      = cos(2*pi*F.*mt-A);
    lambda   =  sum(Amp^2*csq.^2,3)*dfreq*dacro; % non-centrality
    beta     = 1-ncfcdf( finv(1-alpha,2,Nmeas-3) ,2,Nmeas-3,lambda);
    beta_val = sum(beta,"all");
end