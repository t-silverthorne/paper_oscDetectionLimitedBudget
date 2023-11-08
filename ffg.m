N     = 24; % global params
Nfreq = 1e3;
N1    = 4; 
N2    = N-N1;
fmin  = 1;
fmax  = 24;
fvec  = linspace(fmin,fmax,Nfreq)';

tvec1 = linspace(0,1,N1+1)'; % measurement vecs
tvec1 = tvec1(1:end-1);
tvec2 = linspace(0,1,N2+1)';
tvec2 = tvec2(1:end-1);

mloc='seq';
% fminimax(@(tau) -detFvec(getTvec(tau,tvec1,tvec2,mloc),fvec) ...
%                 ,rand,[],[],[],[],0,1)

taus = 0:.001:1;
ys   = arrayfun(@(tau) min(detFvec(getTvec(tau,tvec1,tvec2,mloc),fvec)),taus)
plot(taus,ys)

function tvec = getTvec(tau,tvec1,tvec2,method)
% get parameterization of time vector for reduced design space
switch method
    case 'poly'
        tvec = [mod(tau + tvec1,1); tvec2];
    case 'seq'
        svec1 = tau*tvec1;
        svec2 = (1-tau)*tvec2;
        tvec  = [svec1; tau + svec2];
end
end

function J = detFvec(tvec,fvec)
% evaluate determinant of reduced FIM at all points in fvec
Cmat=cos(2*pi*fvec*tvec');
Smat=sin(2*pi*fvec*tvec');

J=sum(Cmat.*Cmat,2).*sum(Smat.*Smat,2) - sum(Cmat.*Smat,2).^2;
end