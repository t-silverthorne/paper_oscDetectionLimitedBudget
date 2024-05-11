%syms a b phi f ti tj t fmin fmax
%
%% Compute derivatives of NCP
%c2fun   = @(t) cos(2*pi*f*t - phi)^2;
%dc2dphi = @(t) diff(c2fun(t),phi);
%dc2df   = @(t) diff(c2fun(t),f);
%
%% exact expresssion L2norm squraed -- off diagonal terms 
%phi_int = int(dc2dphi(ti)*dc2dphi(tj) + dc2df(ti)*dc2df(tj),phi,0,2*pi);
%f_int = int(phi_int,f,fmin,fmax);
%simplify(f_int)
%
%% exact expresssion L2norm squraed -- diagonal terms 
%phi_int = int(dc2dphi(t)*dc2dphi(t) + dc2df(t)*dc2df(t),phi,0,2*pi);
%f_int = int(phi_int,f,fmin,fmax);
%simplify(f_int)

fmin=1;
fmax=24;

N=randsample(10:100,1);
tau = linspace(0,1,N+1);
tau = tau(1:(end-1));
tau = tau';

Q = ((sin(4*pi*fmax*(tau - tau')) - sin(4*pi*fmin*(tau - tau'))).*(4*tau.*tau'*pi^2 + 1))./(4*(tau - tau'));
Q(eye(size(Q))>0)=pi*(4*pi^2*tau.^2 + 1)*(fmax - fmin);

gersh = NaN(1,N);
for ii=1:N
  gersh(ii) = Q(ii,ii) > sum(abs(Q(ii,[1:(ii-1) (ii+1):N])));
end
prod(gersh)
%f1 = @(ti,tj) ((sin(4*pi*fmax*(ti - tj)) - sin(4*pi*fmin*(ti - tj)))*(4*ti.*tj*pi^2 + 1))./(4*(ti - tj))
%f2 = @(t) pi*(4*pi^2*t.^2 + 1)*(fmax - fmin)
