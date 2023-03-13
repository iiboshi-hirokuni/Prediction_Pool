%% Objective function for hessian computation

function fx = hessn_fcn(para,v)

modelpara = para.*v.pmaskinv + v.pfix.*v.pmask;

[lnpy,retcode,obserror,obsvar] = evaldsge(modelpara,v.data,v.nshock,v.ZZ);

lnprio = priodens(modelpara,v.pmean,v.pstdd,v.pshape);

%fx = lnpy + lnprio;
fx = real(lnpy + lnprio);
