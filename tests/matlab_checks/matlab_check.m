mt  = 0:.1:1;
mt  = mt(1:end-1);
mt  = reshape(mt,length(mt),1)
mt= rand(10,1)
f =2
X =@(delta,f) [ones(length(mt+delta),1) sin(2*pi*f*(mt+delta)) cos(2*pi*f*(mt+delta))]

trF=@(delta,f) trace(X(delta,f)'*X(delta,f))
detF=@(delta,f) det(X(delta,f)'*X(delta,f))


del=0:.1:8;
f=1:.1:8;
[D,F]=meshgrid(del,f);

Z = NaN(size(D));
for ii=1:size(D,1)
    for jj=1:size(D,2)
        Z(ii,jj) = detF(D(ii,jj),F(ii,jj));
    end
end


contourf(D,F,Z,100,'LineStyle','none')