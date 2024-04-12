syms a b phi f ti tj t fmin fmax

% Compute derivatives of NCP
c2fun   = @(t) cos(2*pi*f*t - phi)^2;
dc2dphi = @(t) diff(c2fun(t),phi);
dc2df   = @(t) diff(c2fun(t),f);

% Exact expression for L2 norm with uniform prior
phi_int = int(dc2dphi(ti)*dc2dphi(tj) + dc2df(ti)*dc2df(tj),phi,0,2*pi);
f_int = int(phi_int,f,fmin,fmax);
simplify(f_int)