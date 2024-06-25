clear all
mt   = 0:.1:1;
mt   = mt(1:end-1);
%%
freq = 5.7;
Amp  = 1;
acro = 0;

X     = [ones(length(mt),1) cos(2*pi*freq*mt)' sin(2*pi*freq*mt)'];
Xr    = [cos(2*pi*freq*mt)' sin(2*pi*freq*mt)'];
beta0 = [0 Amp*cos(acro) Amp*sin(acro)];

A=[0 0;
   1 0;
   0 1];

B1 = A'*((X'*X)\A)% formula for B

n   = length(mt);
Ds  = Xr'*Xr;
Bs  = [sum(cos(2*pi*freq*mt)) sum(sin(2*pi*freq*mt))];
rho = n - Bs*(Ds\Bs');
u   = sqrt(rho)*Ds\Bs';
B2 = inv(Ds) + rho*Ds\(Bs'*Bs*inv(Ds)); % Schur complement formula for B

iB1 = inv(B1);
iB2 = Ds - Ds*u*u'*Ds/(1+u'*Ds*u); % Woodbury inverse

max(reshape(abs(B1-B2),4,[]));
max(reshape(abs(iB1-iB2),4,[]));


B3 = inv(Ds-Bs*Bs'/n)

%% Verify NCP formula
clf
tiledlayout(3,1)
freqs=[0.5 1 1.1]
mt   = 0:.02:1;
mt   = mt(1:end-1);
for freq=freqs
    nexttile
    Amp  = 2;
    acro = pi/2;
    Nmc  = 1e4;
    
    X     = [ones(length(mt),1) cos(2*pi*freq*mt)' sin(2*pi*freq*mt)'];
    beta0 = [Amp*cos(acro) Amp*sin(acro)]';
    
    lambda_new = beta0'*((A'*((X'*X)\A) )\beta0);
    lambda_old = sum(Amp^2*cos(2*pi*mt-acro).^2);
    
    
    Ydat = randn([Nmc,length(mt)]) + Amp*cos(2*pi*freq*mt)*cos(acro)+Amp*sin(2*pi*freq*mt)*sin(acro);
    
    betas = (X'*X)\(X'*Ydat');
    
    fits = (X*betas)';
    
    TSS = sum((Ydat-mean(Ydat,2)).^2,2);
    RSS = sum((Ydat-fits).^2,2);
    
    n=length(mt);
    r=3;
    Fstat = (TSS-RSS)*(n-r)/(r-1)./RSS;
    histogram(Fstat,'Normalization','pdf')
    hold on
    
    xv = linspace(0,100,1e4);
    plot(xv,ncfpdf(xv,2,n-3,lambda_new),'-k','LineWidth',3)
    plot(xv,ncfpdf(xv,2,n-3,lambda_old),'-r','LineWidth',3)
    title(strcat('freq = ',string(freq)))
    xlim([0,40])
end
%% verify convexity
%cfun = @(X) min(eig(inv(A'*((X'*X)\A))))
acro = rand*2*pi;
beta = [cos(acro) sin(acro)]';
cfun = @(X) beta'*((A'*((X'*X)\A) )\beta);

Nmeas = randsample(100,1);
mt1 = rand(1,Nmeas);
mt2 = rand(1,Nmeas);
X1 =[ones(length(mt1),1) cos(2*pi*freq*mt1)' sin(2*pi*freq*mt1)'];
X2 =[ones(length(mt2),1) cos(2*pi*freq*mt2)' sin(2*pi*freq*mt2)'];
alpha = rand;
cfun(alpha*X1+(1-alpha)*X2) <= alpha*cfun(X1) + (1-alpha)*cfun(X2)
%%
Nfine=100;

cfun = @(M) min(eig(inv(A'*((M)\A) )))

tau   = (1:Nfine)/Nfine -1/Nfine;
mu1   = randsample([0,1],Nfine,true);
mu2   = randsample([0,1],Nfine,true);
alpha = rand;
LHS   = cfun(getFIM(alpha*mu1+(1-alpha*mu2),tau,freq));
RHS   = alpha*cfun(getFIM(mu1,tau,freq))+(1-alpha)*cfun(getFIM(mu2,tau,freq));
LHS   >= RHS
function M = getFIM(mu,tau,freq)
    mu=reshape(mu,1,1,[]);
    tau=reshape(tau,1,1,[]);
    
    X = [ones(1,1,length(tau)); cos(2*pi*freq*tau); sin(2*pi*freq*tau)];

    M=pagemtimes(X,pagetranspose(X));
    M= sum(M.*mu,3);
end
%%


%%
Ydat = randn(1,length(mt)) + Amp*cos(2*pi*freq*mt-acro);

X = [ones(length(mt),1) cos(2*pi*freq*mt)' sin(2*pi*freq*mt)'];

beta = (X'*X)\(X'*Ydat');
fits = (X*beta)';

TSS = sum((Ydat-mean(Ydat,2)).^2,2);
RSS = sum((Ydat-fits).^2,2);

n=length(mt);
r=3;
Fstat = (TSS-RSS)*(n-r)/(r-1)./RSS
