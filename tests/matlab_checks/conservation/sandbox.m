clear all
clf
mt = rand(1,24+1);%linspace(0,1,24+1);
mt = mt(1:end-1);
mt;

Nfreq = 2^8;
Nacro = 2^8;
fmin  = 1;
fmax  = 24;

freq  = linspace(fmin,fmax,Nfreq);
acro  = linspace(0,2*pi,Nacro+1);
acro  = acro(1:end-1);

[F,A] = meshgrid(freq,acro);

mt       = reshape(mt,1,1,length(mt));
MT       = mt;

Nmeas    = length(mt); % num samples
alpha    = 0.05; % type I error rate
Amp   = 1;


scales=[1 1e-1 1e-2 1 1 1 1 1];
tiledlayout(length(scales),1)
bvec=[];
for sc=scales
    %nexttile
    mt     = reshape(rand(1,20),1,1,20);%MT*sc;
    csq    = cos(2*pi*F.*mt-A);
    lambda =  sum(Amp^2*csq.^2,3); % non-centrality
    beta   = 1-ncfcdf( finv(1-alpha,2,Nmeas-3) ,2,Nmeas-3,lambda);
    %contourf(F,A,beta,200,'LineStyle','none')
    bvec(end+1)=sum(beta,"all");
end
bvec/max(bvec)


%%
ncfSer = @(x,d1,d2,lambda,kk) exp(-lambda/2)*((lambda/2).^kk).*(betainc(d1*x/(d2+d1*x),d1/2+kk,d2))./factorial(kk);
close all
x=rand*2;
lambda=20
Nmeas=10
dlv=.001
lv=0:dlv:10
plot(lv,ncfcdf(x,2,Nmeas-3,lv))
hold on
plot(lv(1:end-1),abs(diff(ncfcdf(x,2,Nmeas-3,lv)))/dlv)
%ncfSer(x,2,Nmeas-3,lambda,0:10)

%%
k=9;
t=rand*10;
s=rand*10;
%integral(@(x) exp(-cos(x-t).^2).*cos(x-t).^(2*k),0,2*pi)

integral(@(x)(cos(x-t).*cos(x-s)).^2,0,2*pi)

%%
syms x t b N c
assume(N>3)
assume(c>0) % check c >0
d1=2;
d2=N-3;
gamma(d1)
betaF=gamma(d1/2)*gamma(d2/2)/gamma(d1/2+d2/2);
limit(simplify(int((1-t)^(d2/2-1),0,d1*x/(d2+d1*x)))/betaF,N,Inf)

%%
syms x
int(cos(x)^2*sin(x)^2,x,0,2*pi)