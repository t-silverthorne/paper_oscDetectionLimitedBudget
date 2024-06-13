%% is phase average conserved
clear all
alpha=.05;
Nacro=2^6;
Nmeas=20;

clf
%%
freq=10;
mt=rand(1,Nmeas);
acro  = linspace(0,2*pi,Nacro+1);
acro  = acro(1:end-1);
dacro = 2*pi/Nacro;

Amp=1.1;


%mt = linspace(0,1,Nmeas+1);
%mt = [0.1*mt(1:end-1) .5];

acro = reshape(acro,1,[]);
mt   = reshape(mt,[],1);
csq      = cos(2*pi*freq*mt-acro);
lambda   = sum(Amp^2*csq.^2,1);
beta     = 1-ncfcdf( finv(1-alpha,2,Nmeas-3) ,2,Nmeas-3,lambda);
mean(beta)
plot(acro,beta)
hold on