m = 16; n = 8;
A = randn(m,n);
b = randn(m,1);
p = 4;
C = randn(p,n);
d = randn(p,1);

cvx_begin
    variable x(n) integer;
    minimize( norm(A*x-b) );
    subject to
        C*x <= d;
        norm(x,Inf) <= 5;
cvx_end

%% similar to E-optimality
clear
n = 10
V = randn(n,3)
cvx_begin
    variable x(n) binary;
    minimize(trace(V'*diag(x)*V))
    subject to
        sum(x) == 5
cvx_end
x
cvx_optval
%%
clear
n = 10;
%V = randn(n,2);
V = [1:2:20; 2:2:20]'


A=[0 0;
   1 0;
   0 1];
cvx_begin
    variable x(n) binary;
    minimize(lambda_max(V'*diag(x)*V))
    subject to
        sum(x) == 5
cvx_end
