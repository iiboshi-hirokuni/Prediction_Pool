function [rloglh,retcode,obsmean,obsvar] = ...
    evaldsge_dr(para, yy, sysmat, HH, psi, zz1, casefile);

npara   = size(para, 1);
retcode = 0;

% solve the DSGE model

[T1,TC,TEPS,RC] = dsgesolv(para, sysmat);

nseries  = size(yy, 2);
nstate   = size(T1, 2)*2;
nshock = size(TEPS, 2);

nobs      = size(yy, 1);
nvar      = size(yy, 2);
loglh     = 0;
loglhzero = -1E8;
obsmean   = zeros(nobs, nseries);
obsvar    = zeros(nobs, nseries);

% -------------------
% Check determinacy 
%--------------------
if (RC(1) == 1) && (RC(2)==1);
   %/* determinacy */
   retcode(1) = 0;
   TT = T1;
   RR = TEPS;
   
elseif (RC(1) == 1) && (RC(2)==0) 
   %/* indeterminacy */
   retcode(1) = 1;
   TT = T1;
   RR = TEPS;
   rloglh = loglhzero;
   return;

else
   %/* no equilibrium exists, numerical problems */
   retcode(1) = RC(1);
   rloglh = loglhzero;
   return;

end

% create system matrices for state space model

% These matrices are regime independent
%zz = make_zz(casefile, pv.nstate); 
%gam = zz;

DD = zeros(nvar, 1);      
QQ = createcov(para(npara-nshock+1:npara,1));
VV = zeros(nshock, nvar);

TT = [T1,zeros(nstate/2,nstate/2);...
     eye(nstate/2),zeros(nstate/2,nstate/2)];

RR = [TEPS;...
      zeros(nstate/2,9)];
ZZ = [zz1, (-1)*diag(psi)*zz1];

% Check whether covariance matrix QQ is positive definite

if sum(eig(QQ) < 0) > 0
   loglh = loglhzero;
   return;
end

% We can now define the initial mean and variance for the state vector

At = make_init(yy, casefile, nstate);
Pt = dlyap(TT,RR*QQ*RR');

% compute likelihood with Kalman filter

t = 1;
while t <= nobs
   
   At1 = At;
   Pt1 = Pt;
   
   % Forecasting
   alphahat = TT*At1;
   Phat = TT*Pt1*TT' + RR*QQ*RR';
   yhat = ZZ*alphahat + DD;
   nu   = yy(t,:) - yhat';
   
   Ft   = ZZ*Phat*ZZ' + HH + ZZ*RR*VV + (ZZ*RR*VV)';
   Ft   = 0.5*(Ft + Ft');
   
   loglh = loglh -0.5*size(yy, 2)*log(2*pi)-0.5*log(det(Ft)) ...
           - 0.5*nu*inv(Ft)*nu';
   
   % Updating
   At = alphahat + (Phat*ZZ' + RR*VV)*inv(Ft)*nu';
   Pt = Phat - (Phat*ZZ'+RR*VV)*inv(Ft)*(Phat*ZZ'+RR*VV)';
   
   %  store
   obsmean(t,:) = yhat';
   obsvar(t,:)  = diag(Ft)';
   
   t = t+1;
end  

rloglh = real(loglh);

%-----------------------------
function omega = createcov(para)

npara = max(size(para));
omega = zeros(npara, npara);
for i = 1:npara
  omega(i, i) = para(i)^2;
end



