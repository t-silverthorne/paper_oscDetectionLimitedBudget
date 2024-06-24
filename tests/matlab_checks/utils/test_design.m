% mt=[0.000000000 0.006944444 0.020833333 0.034722222 0.041666667 0.048611111 ...
% 0.055555556 0.062500000 0.069444444 0.076388889 0.083333333 0.090277778 ...
% 0.097222222 0.104166667 0.111111111 0.118055556 0.125000000 0.131944444 ...
% 0.138888889 0.145833333 0.152777778 0.159722222 0.166666667 0.180555556 ...
% 0.187500000 0.194444444 0.201388889 0.215277778 0.236111111 0.250000000 ...
% 0.270833333 0.284722222];
Nmeas=32;
mt=linspace(0,.25,Nmeas+1);
mt=mt(1:end-1);
Nmeas=length(mt);
Amp=0.85;
acro=0;
freq=3.95;

nrep=1e2

freqs = 3.5:.1:4;
cmap=cool(length(freqs));
Nmc = 1e2;
clf
for ii=1:length(freqs)
    freq=freqs(ii);
    pwr_exact = eval_exact_power(Amp,acro,freq,mt)
    
    nreps = [1e1 1e2 1e3 1e4 1e5];
    err = [];
    for nrep=nreps
        
        
        pval=[];
        parfor jj=1:nrep
            Y=randn(Nmc,Nmeas) + Amp*cos(2*pi*freq*mt-acro);
            res=fit_cosinor_model(Y,mt,1/freq,'householder');
            pval = [pval; res.pval];
        end
        pwr_est = mean(pval < .05);
        err(end+1)= abs(pwr_est-pwr_exact)/pwr_exact
    end
    loglog(Nmc*nreps,err,'Color',cmap(ii,:))
    hold on
    drawnow
    pause(0.1)
end
