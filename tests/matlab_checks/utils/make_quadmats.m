function  Amlist = make_quadmats(opts)
N=opts.N;Nfreq=opts.Nfreq;fmin=opts.fmin;fmax=opts.fmax;% unpack

tvec=(0:N)/N;% make tvec
tvec=tvec(1:end-1);
length(tvec);

fvec = linspace(fmin,fmax,Nfreq); % make fvec
fvec = reshape(fvec,1,1,length(fvec));
tvec = reshape(tvec,length(tvec),1);

cvec = cos(2*pi*tvec.*fvec);% to tensor
svec = sin(2*pi*tvec.*fvec);

a11 = cvec.*cvec;% intermediate matrices
a12 = cvec.*svec;
a22 = svec.*svec;
a11T = pagetranspose(a11);
a12T = pagetranspose(a12);
a22T = pagetranspose(a22);

% result
Amlist = pagemtimes(a11,a11T) + pagemtimes(a22,a22T) - ...
    ( pagemtimes(a11,a22T) + pagemtimes(a22,a11T)) + ...
    4*pagemtimes(a12,a12T);
end