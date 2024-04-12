N = 1e2;
u = rand(N,1);
v = randn(N,1);
w = rand(N,1);

x = rand(N,1);

%%
A= [x'*u x'*v;
    x'*v x'*w];

B = u*u' + w*w' - (u*w' + w*u') + 4*v*v';

lmin_1 = min(eig(A));


lmin_2  = 0.5*(trace(A) - sqrt( x'*B*x));

%%
eig(B)