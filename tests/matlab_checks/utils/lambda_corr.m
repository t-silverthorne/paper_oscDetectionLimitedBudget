clear
global freq Amp acro
freq=.65;
Amp=2;
acro=pi/2;

mt   = 0:.1:1;
mt   = mt(1:end-1);
lambda=sum(Amp^2*cos(2*pi*freq*mt-acro).^2);
n=length(mt);
Nmc=1e4;
clf
histogram(arrayfun(@(ind)get_num(mt),1:Nmc),'Normalization','pdf')
xv = linspace(0,100,1e3);
hold on
plot(xv,ncx2pdf(xv,2,lambda),'-k')

%%
freq=1.7;
Amp=1;
acro=pi/2;
Ydat = randn(1,length(mt)) + Amp*cos(2*pi*freq*mt-acro);

X = [ones(length(mt),1) cos(2*pi*freq*mt)' sin(2*pi*freq*mt)'];

beta = (X'*X)\(X'*Ydat');
fits = (X*beta)';

TSS = sum((Ydat-mean(Ydat,2)).^2,2);
RSS = sum((Ydat-fits).^2,2);

A=[0 0;
   1 0;
   0 1];
A'*beta;

B=A'*(pinv(X'*X)*A)

Xr=[cos(2*pi*freq*mt)' sin(2*pi*freq*mt)'];

Xr'*Xr
inv(B)

%%

freq=.4;
Amp=1;
acro=pi/2;
Ydat = randn(1,length(mt)) + Amp*cos(2*pi*freq*mt-acro);

X = [ones(length(mt),1) cos(2*pi*freq*mt)' sin(2*pi*freq*mt)'];

beta = (X'*X)\(X'*Ydat');
fits = (X*beta)';

TSS = sum((Ydat-mean(Ydat,2)).^2,2);
RSS = sum((Ydat-fits).^2,2);

A=[0 0;
   1 0;
   0 1];


B=A'*(pinv(X'*X)*A);

TSS-RSS
num=(A'*beta)'*(B\(A'*beta))

lamba_naive = sum(Amp^2*cos(2*pi*freq*mt-acro).^2)

beta0=[0 Amp*cos(2*pi*acro) Amp*sin(2*pi*acro)]';
B=A'*(pinv(X'*X)*A);
lambda = (A'*beta0)'*(B\(A'*beta0));

%% try ncfpdf with new lambda

Nmc=5e5

freq = 0.6;
Amp  = 2;
acro = pi/5;
Ydat = randn([Nmc,length(mt)]) + Amp*cos(2*pi*freq*mt)*cos(acro)+Amp*sin(2*pi*freq*mt)*sin(acro);

A=[0 0;
   1 0;
   0 1];
X = [ones(length(mt),1) cos(2*pi*freq*mt)' sin(2*pi*freq*mt)'];
beta0=[0 Amp*cos(acro) Amp*sin(acro)]';
B=A'*((X'*X)\A);

lambda_naive = sum(Amp^2*cos(2*pi*freq*mt-acro).^2)
lambda       = (A'*beta0)'*(B\(A'*beta0))

trace(B\A'*beta0*beta0'*A)

%U=chol(B);
%lambda_alt = sum((U'\(A'*beta0)).^2)

betas = (X'*X)\(X'*Ydat');

fits = (X*betas)';

TSS = sum((Ydat-mean(Ydat,2)).^2,2);
RSS = sum((Ydat-fits).^2,2);

n=length(mt);
r=3;
Fstat = (TSS-RSS)*(n-r)/(r-1)./RSS;
clf
histogram(Fstat,'Normalization','pdf')
hold on

plot(xv,ncfpdf(xv,2,n-3,lambda),'-k','LineWidth',3)
plot(xv,ncfpdf(xv,2,n-3,lambda_naive),'-r','LineWidth',3)
xlim([0,40])
%%
B=[sum(cos(2*pi*freq*mt)); sum(sin(2*pi*freq*mt))]
C=B'
A=length(mt)
D= [cos(2*pi*freq*mt)*cos(2*pi*freq*mt)' cos(2*pi*freq*mt)*sin(2*pi*freq*mt)';
    cos(2*pi*freq*mt)*sin(2*pi*freq*mt)' sin(2*pi*freq*mt)*sin(2*pi*freq*mt)']
D
%%
A-B'*(D\C')

%%
clf
dat=A'*betas;
histogram2(dat(1,:)',dat(2,:)','Normalization','pdf')
hold on

x=-3:.01:3
y=-3:.01:3
[X,Y] = meshgrid(x,y);
Z=mvnpdf([X(:) Y(:)],(A'*beta0)',B);
Z=reshape(Z,length(y),length(x));
surf(X,Y,Z,'LineStyle','none','FaceAlpha',.75)



%%
function num=get_num(mt)
global freq Amp acro

Ydat = randn(1,length(mt)) + Amp*cos(2*pi*freq*mt-acro);


X = [ones(length(mt),1) cos(2*pi*freq*mt)' sin(2*pi*freq*mt)'];

beta = (X'*X)\(X'*Ydat');
fits = (X*beta)';

TSS = sum((Ydat-mean(Ydat,2)).^2,2);
RSS = sum((Ydat-fits).^2,2);

A=[0 0;
   1 0;
   0 1];
A'*beta;

B=A'*(pinv(X'*X)*A);

TSS-RSS
num=(A'*beta)'*(B\(A'*beta))

end



