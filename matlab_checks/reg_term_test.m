syms a b phi f ti tj t fmin fmax

% Compute derivatives of NCP
c2fun   = @(t) cos(2*pi*f*t - phi)^2;
dc2dphi = @(t) diff(c2fun(t),phi);
dc2df   = @(t) diff(c2fun(t),f);

% exact expresssion L2norm squraed -- off diagonal terms 
phi_int = int(dc2dphi(ti)*dc2dphi(tj) + dc2df(ti)*dc2df(tj),phi,0,2*pi);
f_int = int(phi_int,f,fmin,fmax);
simplify(f_int)

% exact expresssion L2norm squraed -- on diagonal terms 
phi_int = int(dc2dphi(t)*dc2dphi(t) + dc2df(t)*dc2df(t),phi,0,2*pi);
f_int = int(phi_int,f,fmin,fmax);
simplify(f_int)