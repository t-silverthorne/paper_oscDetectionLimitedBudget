N  = datasample(3:10,1);
mt = linspace(0,1,N);
X  = [ones(N,1) 50*randn(N,2)];
%X  = [ones(N,1) 2*pi*mt' 2*pi*mt'.^2];
%X  = [ones(N,1) cos(2*pi*mt') sin(2*pi*mt')];
r  = rank(X)
U=eye(N)-X*pinv(X);
V=(eye(N)-ones(N,N)/N);

V-U
(V-U)^2

rank(V-U,1e-12)
r-1
%%
close all
N     = 10;
Nsamp = 1e4;
ra    = datasample(5:N,1);
rb    = datasample(1:(ra-3),1);
A     = rand(N,N);
B     = rand(N,N);

[Qa,~]=qr(A);
[Qb,~]=qr(B);

U=Qa(1:ra,:)'*Qa(1:ra,:);
V=Qb(1:rb,:)'*Qb(1:rb,:);

max(max(abs(U^2-U)))
max(max(abs(V^2-V)))

x=randn(N,Nsamp);
xdat=diag(x'*(U-V)*x);

histogram(xdat, 'Normalization', 'pdf');
hold on;
ax = gca;
x = linspace(ax.XLim(1), ax.XLim(2), 1000);
plot(x, chi2pdf(x,ra-rb), 'LineWidth', 2);
hold off;

%%

% rank(U,1e-12)
% N-r

y=rand(N,1);

TSS = norm(y-mean(y)*ones(N,1),2)^2
RSS = norm(( eye(N) - X*((X'*X)\X') )*y,2)^2

TSS-RSS
y'*(V'*V-U'*U)*y

%%
norm((eye(N)-ones(N,N)/N - U)*y,2)
%%

V=eye(N)-ones(N,N)/N;
rank(V-U)


