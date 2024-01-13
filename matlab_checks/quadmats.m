addpath('utils/')
%% Linfty test
N     = 5;
Nfreq = 3;
fmin  = 1;
fmax  = 24;
Amlist = make_quadmats(opts);

Amlist(:,:,1)
Amlist(:,:,2)
Amlist(:,:,3)

%% L1 test
opts=struct();
opts.N     = 8;
opts.Nfreq = 64;
opts.fmin  = 1;
opts.fmax  = 24;
Amlist = make_quadmats(opts);

mean(Amlist,3)
