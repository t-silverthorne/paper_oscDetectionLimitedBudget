syms f t phi ti tj fmin fmax
assume(f,'real')
assume(phi,'real')
assume([t ti tj],'real')
ncp = cos(2*pi*f*t-phi)^2;

dncp_df   = @(s) subs(diff(ncp,f),t,s)
dncp_dphi = @(s) subs(diff(ncp,phi),t,s)
int(dncp_df(ti)*dncp_df(tj) + dncp_dphi(ti)*dncp_dphi(tj),f,fmin,fmax)
%%
int(int(dncp_df(ti)*dncp_df(tj),phi,0,2*pi),f,fmin,fmax)

int(int(dncp_dphi(ti)*dncp_dphi(tj),phi,0,2*pi),f,fmin,fmax)

%%
int(gr'*gr,phi,0,2*pi)
%%
phi_num=0
fu = @(f,t) cos(2*pi*f*t-phi_num)^2

f0 = 1.03%10*rand
h  = 1e-8
t0 = 0.005%10*rand

(fu(f0*(1+h),t0)-fu(f0,t0))/f0/h

eval(subs(diff(ncp,f),[f,phi,t],[f0,phi_num,t0]))