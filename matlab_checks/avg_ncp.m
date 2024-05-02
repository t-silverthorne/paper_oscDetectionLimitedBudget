syms a b phi f ti tj t fmin fmax

%%
int(cos(2*pi*f*t-phi)^2,phi,0,2*pi)
%%
simplify(int(int(cos(2*pi*f*t-phi)^2,phi,0,2*pi),f,fmin,fmax)/2/pi/(fmax-fmin))
%%
solve(diff(simplify(int(cos(2*pi*f*t-phi)^2,f,fmin,fmax)),phi),phi)
