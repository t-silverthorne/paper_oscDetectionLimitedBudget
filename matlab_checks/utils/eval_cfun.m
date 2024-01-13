function fval = eval_cfun(x,Amat,opts)
fval=NaN;
if strcmp(opts.costfun_type,'L1')
    fval=x'*Amat*x;
elseif strcmp(opts.costfun_type,'Linfty')
    fval=max(pagemtimes(x',pagemtimes(Amat,x)));
end
end

