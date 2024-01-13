addpath('utils/')
opts=struct();
opts.N     = 2^4;
opts.Nfreq = 2^6;
opts.fmin  = 1;
opts.fmax  = 24;

Amlist = make_quadmats(opts);
x      = [0 1 0 0 0 0 1 0 1 0 0 0 0 0 1 0];
x      = reshape(x,length(x),1);

opts.costfun_type='L1';
eval_cfun(x,mean(Amlist,3),opts)

opts.costfun_type='Linfty';
eval_cfun(x,Amlist,opts)



