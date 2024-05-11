function worstpwr = eval_worst_power(mt,freqs,acros,Amp,alpha)
% assumes that mt, freqs and acros already are of the following shapes
%   mt     = reshape(mt,[length(mt) 1 1 1]);
%   acros  = reshape(acros,[1 1 length(acros) 1]);
%   freqs  = reshape(freqs,[1 1 1 length(freqs)]);
Nmeas  = length(mt);
Xarg   = 2*pi*freqs.*mt-acros;
csq    = cos(Xarg);
lambda = Amp^2*pagemtimes(pagetranspose(csq),csq);
beta   = 1-ncfcdf( finv(1-alpha,2,Nmeas-3) ,2,Nmeas-3,lambda);

worstpwr = min(beta);
end




