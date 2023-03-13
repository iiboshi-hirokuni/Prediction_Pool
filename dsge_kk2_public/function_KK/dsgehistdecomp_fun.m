function [yy_decomp,yy_smooth1,yy_smooth2,DD,leverage] = dsgehistdecomp_fun(para,nforecast,nvar,nshock);
%/* INPUT:
%**   para   : parameter vector 
%**   nforecast   : length of forecast
%** OUTPUT:
%**   ssforecast  : impulse response 
%**            the structural shocks.
%**            0 if retcode=0
%**            (nforecast by 25)
%*/

global ZZ YY ;


% solve DSGE

[T1,TC,T0,RC] = dsgesolv(para);


 [yy_decomp,yy_smooth1,yy_smooth2,DD,leverage] = Smooth_dsge(para,YY,nshock,ZZ,T1,TC,T0);




