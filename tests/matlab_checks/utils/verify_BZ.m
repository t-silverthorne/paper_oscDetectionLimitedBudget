sim_choice='infradian';
Nmc  = 1e4;
acro = pi/2;
Amp  = 0;
switch sim_choice
    case 'infradian'
        mt   = 0:.1:1;
        mt   = mt(1:end-1);
        freq = 0.5;
    case 'circadian'
        mt   = 0:.1:1;
        mt   = mt(1:end-1);
        freq = 1.1;
end

%% verify Lemma 1
Amp = 0
Ydat = randn(Nmc,length(mt)) + Amp*cos(2*pi*freq*mt-acro);

X = [ones(length(mt),1) cos(2*pi*freq*mt)' sin(2*pi*freq*mt)'];

betas = (X'*X)\(X'*Ydat');
fits = (X*betas)';

TSS = sum((Ydat-mean(Ydat,2)).^2,2);
RSS = sum((Ydat-fits).^2,2);

n=length(mt);
r=3;
Fstat = (TSS-RSS)*(n-r)/(r-1)./RSS;

xv = linspace(0,100,1e3);
clf
histogram(Fstat,'Normalization','pdf')
hold on
plot(xv,fpdf(xv,2,n-3),'-k')
xlim([0 10])

%% verify Lemma 2
Nmc=1e6
A=eye(2);
freq=0.5
Amp=1;
acro=pi/2
Ydat = randn(Nmc,length(mt)) + Amp*cos(2*pi*freq*mt-acro);

X = [ones(length(mt),1) cos(2*pi*freq*mt)' sin(2*pi*freq*mt)'];

betas = (X'*X)\(X'*Ydat');

fits = (X*betas)';

TSS = sum((Ydat-mean(Ydat,2)).^2,2);
RSS = sum((Ydat-fits).^2,2);



betas=betas(2:3,:);
X = [cos(2*pi*freq*mt)' sin(2*pi*freq*mt)'];
eig(pinv(X'*X))
%%
betas=reshape(betas,2,1,[]);

lambda = sum(Amp^2*cos(2*pi*freq*mt-acro).^2)
beta_true = [Amp*cos(acro); Amp*sin(acro)];
beta_true'*(pinv(X'*X)\beta_true)

% Individual parts of lemma 3

nums=reshape(pagemtimes(pagetranspose(betas),pagemldivide(pinv(X'*X),betas)),Nmc,1,[]);
clf
tiledlayout(2,1)
nexttile(1)
histogram(nums,'Normalization','pdf')
hold on
plot(xv,ncx2pdf(xv,2,lambda))
ylim([0 0.1])

nexttile(2)
histogram((TSS-RSS),'Normalization','pdf')
hold on
plot(xv,ncx2pdf(xv,2,lambda))
ylim([0 0.1])

%% Full Lemma 3
lambda_est = reshape(pagemtimes(pagetranspose(betas),pagemldivide(pinv(X'*X),betas)),...
    Nmc,1,[]);
clf
histogram(lambda_est*(n-r)./RSS/(r-1),'Normalization','pdf')
%histogram((TSS-RSS)*(n-r)/(r-1)./RSS,'Normalization','pdf')
hold on
plot(xv,ncfpdf(xv,2,n-3,lambda))
xlim([0 100])
ylim([0 0.1])
