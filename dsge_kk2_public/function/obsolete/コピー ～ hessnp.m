% Compute Hessian

function [hessian,errorij] = hessn(fcn,para,v)

npara = size(para,1);
ndx = 6;
for i = 1:ndx
    dx(i) = 5 + 2*(i-1);
end
dx = exp(-dx);
hessian = zeros(npara,npara);
%gradx    = zeros(ndx,1);
%grady    = zeros(ndx,1);
%gradxy   = zeros(ndx,1);
%hessdiag = zeros(ndx,1);
dxscale  = ones(npara,1);

fx  =feval(fcn,para,v);

fmt = '%5.0f ';

matlabpool open 2
disp(' ');
disp('Computing Hessian..');
disp('--------------------------------------------------');
disp(sprintf(['  diagonal elements     :   ' fmt ],npara));
disp(sprintf(['  off-diagonal elements :   ' fmt ],npara*(npara-1)/2));

for seli = 1:(npara-1)
   parfor selj =  seli:npara;
      hessdiag = fhessdiag(seli, selj, fcn, para, fx, dx, ndx, dxscale, v);
      hessian(seli,selj) = -0.5*(hessdiag(3)+hessdiag(4));
      if mod(seli+(selj-1)*npara,10) == 0
         disp(sprintf(['                            ' fmt fmt],seli,selj));
      end
   end
end

hessdiag = fhessdiag(npara, npara, fcn, para, fx, dx, ndx, dxscale, v);
hessian(npara,npara) = -0.5*(hessdiag(3)+hessdiag(4));


matlabpool close

errorij = [];

for seli = 1:npara-1
   for selj =  seli+1:npara;
      if ( hessian(seli,seli) == 0 ) || (hessian(selj,selj) == 0)
        corrij = 0;
      else
        corrij = hessian(seli,selj)/sqrt(hessian(seli,seli)*hessian(selj,selj));
      end
      
      if (corrij < -0.98) || (corrij > 0.98);
         hessian(seli,selj)=0.9*sqrt(hessian(seli,seli)*hessian(selj,selj));
         errorij = [errorij ;[seli,selj,corrij]];
      elseif (corrij > -0.005) && (corrij < 0.005);
        hessian(seli,selj)=0;
      end 
      hessian(selj,seli) = hessian(seli,selj);
   end
end

hessian = real(hessian);

function ret = fhessdiag(seli, selj, fcn, para, fx, dx, ndx, dxscale, v)
    ret = zeros(ndx,1);
    for i = 2:3
        paradx = para;
        parady = para;
        paradx(seli) = paradx(seli) + dx(i)*dxscale(seli);
        parady(selj) = parady(selj) - dx(i)*dxscale(selj);
        paradxdy = paradx;
        paradxdy(selj) = paradxdy(selj) - dx(i)*dxscale(selj);
        fdx = feval(fcn,paradx,v);
        fdy = feval(fcn,parady,v);
        fdxdy = feval(fcn,paradxdy,v);
        ret(i) = -( fx - fdx - fdy + fdxdy )/(dx(i)*dx(i)*dxscale(seli)*dxscale(selj));
    end

