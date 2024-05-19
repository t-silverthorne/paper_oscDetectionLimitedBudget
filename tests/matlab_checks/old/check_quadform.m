% setup
N    = 5;
f    = 10*rand;
t    = sort(rand(N,1));
cvec = cos(2*pi*f*t);
svec = sin(2*pi*f*t);

% hadamard vectors
a11 = cvec.*cvec;
a12 = cvec.*svec;
a22 = svec.*svec;

%% Two symmetric matrices
A1  = 0.5*(a11*a22' + a22*a11') - a12*(a12');
A2  = a11*a11' + a22*a22' - a11*a22' - a22*a11' + 4*a12*(a12');

%% check they give same result
mu = randsample([0,1],N,true)';
sum(mu)^2- 4*mu'*A1*mu ==mu'*A2*mu

%% compare spectra
max(abs(eig(A1)-eig(A2)))


%min(eig([cvec'*cvec cvec'*svec; cvec'*svec svec'*svec])) % min eig FIM
%%x'*(u*v' + v*u')*x 
%x'*(2*u*v')*x 