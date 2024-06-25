clear all
mt   = 0:.05:1;
mt   = mt(1:end-1);

freq = 0.7;
Amp  = 1;
acro = 0;

X     = [ones(length(mt),1) cos(2*pi*freq*mt)' sin(2*pi*freq*mt)'];
Xr    = [cos(2*pi*freq*mt)' sin(2*pi*freq*mt)'];
beta0 = [0 Amp*cos(acro) Amp*sin(acro)];

A=[0 0;
   1 0;
   0 1];
min(eig(A'*((X'*X)\A)))

D = Xr'*Xr
%%
b = [sum(cos(2*pi*freq*mt)); sum(sin(2*pi*freq*mt))];
n = length(mt);
1/max(eig(D-b*b'/n))

%%
D
%b*b'/n