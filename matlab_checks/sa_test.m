addpath('utils')

alpha = .05;
Nmeas = 16;
fmin  = 1;
fmax  = 24;
Nf    = 2^10;
Nacro = 2^8;
Amp   = 1.5;
mt    = (0:Nmeas)/Nmeas;
mt    = mt(1:end-1);
freqs = linspace(fmin,fmax,Nf);
acros = linspace(0,2*pi,Nacro+1);
acros = acros(1:end-1);
mt    = reshape(mt,[length(mt) 1 1 1]);
acros = reshape(acros,[1 1 length(acros) 1]);
freqs = reshape(freqs,[1 1 1 length(freqs)]);


eval_worst_power(mt,freqs,acros,Amp,alpha)


sa_cfun = @(x) 1-eval_worst_power(x,freqs,acros,Amp,alpha);


options = optimoptions(@simulannealbnd,'MaxIterations',1000);

num_ens = 1e2;
for ii=1:num_ens
    x0=rand(Nmeas,1);
    x=simulannealbnd(sa_cfun,x0,zeros(Nmeas,1),ones(Nmeas,1),options);
    Amat(ii,:)=x'
end

save('../results/matlab_opt')
%xx=[0.002 0.031 0.180 0.266 0.339 0.415 0.430 0.476 0.483 0.526 0.538 0.598 0.622 0.636 0.905 0.978]
%xx=[0:10]/10;
sa_cfun(x)