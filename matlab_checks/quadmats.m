%% Linfty test
N     = 5;
Nfreq = 3;
fmin  = 1;
fmax  = 24;
Amlist = getAmlist(N,Nfreq,fmin,fmax);

Amlist(:,:,1)
Amlist(:,:,2)
Amlist(:,:,3)

%% L1 test
N     = 8;
Nfreq = 64;
fmin  = 1;
fmax  = 24;
Amlist = getAmlist(N,Nfreq,fmin,fmax);

sum(Amlist,3)

function  Amlist = getAmlist(N,Nfreq,fmin,fmax)
% make tvec
tvec=(0:N)/N;
tvec=tvec(1:end-1);
length(tvec);

% make fvec
fvec = linspace(fmin,fmax,Nfreq);
fvec = reshape(fvec,1,1,length(fvec));
tvec = reshape(tvec,length(tvec),1);

% to tensor
cvec = cos(2*pi*tvec.*fvec);
svec = sin(2*pi*tvec.*fvec);

% intermediate matrices
a11 = cvec.*cvec;
a12 = cvec.*svec;
a22 = svec.*svec;

% page transpose
a11T = pagetranspose(a11);
a12T = pagetranspose(a12);
a22T = pagetranspose(a22);

% result
Amlist = pagemtimes(a11,a11T) + pagemtimes(a22,a22T) - ...
    ( pagemtimes(a11,a22T) + pagemtimes(a22,a11T)) + ...
    4*pagemtimes(a12,a12T);
end