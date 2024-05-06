%% what happens if measuring at single point in time
clf
t=0:.01:2;
f=1:.01:3;

[T,F]=meshgrid(t,f);

ncp = @(t,f) cos(2*pi*t.*f).^2;

g = @(t,f,delta) ncp(t+delta,f) - ncp(t,f);

N = ncp(T,F);

tiledlayout(1,2)
nexttile(1)
contourf(T,F,N,100,'LineStyle','none')

nexttile(2)
dncp = @(t,f) -4*pi*t.*cos(2*pi*f.*t).*sin(2*pi*f.*t);
dg = @(t,f,delta) dncp(t+delta,f)-ncp(t,f);
dN = dncp(T,F)
contourf(T,F,dN,100,'LineStyle','none')
%%
clf
plot(t,g(t,2.01,1))

%% what happens if measuring at several points in time
clf

mt  = 0:.1:1;
mt  = mt(1:end-1);
mt  = reshape(mt,1,1,length(mt));

ncp = @(delta,f) sum(cos(2*pi*(mt+delta).*f).^2,3)
plot(ncp(0:.01:1,1))


del=0:.01:8;
f=1:.01:8;
[D,F]=meshgrid(del,f);

Z=ncp(D,F)
contourf(D,F,Z,100,'LineStyle','none')
%%
clf
plot(del,ncp(del,5))
hold on
plot(del,ncp(del,5.2))
%% several points in time mean over f
mt=mt
fmin =1
fmax =2
h =  @(delta,t) sum((fmax-fmin)/2 + sinc(4*pi*fmax*(t+delta)/pi)*fmax/2- sinc(4*pi*fmin*(t+delta)/pi)*fmin/2,3)
plot(del,h(del,mt))

[u,v]=findpeaks(h(del,mt))
mod(del(v),1)
%% checing
h(0,0)

fmax=1.5
fmin=2.5
h1 = @(t) ((fmax*t)/2 - (fmin*t)/2)/t + (sin(4*pi*fmax*t)/8 - sin(4*pi*fmin*t)/8)/(t*pi)
h2 = @(t) (fmax-fmin)/2 + sinc(4*pi*fmax*t/pi)*fmax/2- sinc(4*pi*fmin*t/pi)*fmin/2

h2(0.1)