function [st_save, smooth1, smooth2 , DD,leverage] = Smooth_dsge(para,yy,nshock,ZZ,T1,TC,TEPS)
% function [st_save, obsmean, DD] = Smooth_dsge(para,yy,nshock,ZZ,T1,TC,TEPS)

global ZZ_leverage;

% smooth_type = 'state_smooth';
smooth_type = 'disturbance_smooth' ;

npara   = size(para, 1);
retcode = 0;

% solve the DSGE model
% 
[T1,TC,TEPS,RC] = dsgesolv(para);
%DD = zeros(size(ZZ, 1),1);

 TT = T1;
 RR = TEPS;

nseries  = size(yy, 2);
nstate   = size(T1, 2);

nobs      = size(yy, 1);
loglh     = 0;
loglhzero = -1E8;
% obsmean   = zeros(nobs, 2*nstate);
% obssmooth =zeros(nobs, 2*nstate);
obsmean =zeros(nobs, nseries);
smooth1 =zeros(nobs, nseries);
smooth2 =zeros(nobs, nseries);
leverage =zeros(nobs,1);
obsvar    = zeros(nobs, nseries);
nu_save   = zeros(nseries,nobs);
ft_save   = zeros(nseries,nseries*nobs); 
kg_save   = zeros(nstate,nseries*nobs);
at_save   = zeros(nobs,2*nstate);
pt_save   = zeros(2*nstate,2*nstate*nobs);
shock     = zeros(nobs, nshock);
st_save = zeros(nobs,nseries*nshock);

% create system matrices for state space model

% These matrices are regime independent


z_star_bar = para(14,1);
psi_bar    = para(15,1);

l_bar   = mean(yy(:,5));  %1.0;
r_n_bar = mean(yy(:,8));
pi_bar  = 0.25;

DD = zeros(nseries,1);
% DD(1, 1) = z_star_bar;
% DD(2, 1) = z_star_bar;
% DD(3, 1) = z_star_bar+psi_bar;
% DD(4, 1) = z_star_bar;
% DD(5, 1) = l_bar;
% DD(6, 1) = pi_bar;
% DD(7, 1) = -psi_bar;
% DD(8, 1) = r_n_bar;
% DD(9, 1) = 0;
% DD(10, 1) = z_star_bar;

HH = zeros(nseries,nseries);
QQ = createcov(para(31:41,1));
VV = zeros(nshock,nseries);

% Check whether covariance matrix QQ is positive definite

if sum(eig(QQ) < 0) > 0
   loglh = loglhzero;
   return;
end

% We can now define the initial mean and variance for the state vector
%
%At = zeros(nstate,1);
At = zeros(nstate*2,1);
% At(14)=yy(1,1); At(2)=yy(1,2); At(15)=yy(1,3); At(3)=yy(1,4);    
% At(30)=yy(1,5); At(11)=yy(1,6); % At(23)= yy(1,7); 
% At(17)=yy(1,8); %At(14)=yy(1,9); 
% At(4)= yy(1,10); % At(14+40)=yy(1,1); At(2+40)=yy(1,2); At(15+40)=yy(1,3); At(3+40)=yy(1,4);    
% At(30+40)=yy(1,5); At(11+40)=yy(1,6); %At(23+40)= yy(1,7); 
% At(17+40)=yy(1,8); %At(14+40)=yy(1,9); At(4+40)= yy(1,10); 

a1 = At;

Pt = dlyap(TT,RR*QQ*RR');
Pt = [Pt zeros(nstate);...
      zeros(nstate) Pt];

% Pt = 0.0001*Pt;
  
p1 = Pt;   

TT = [TT zeros(nstate);
      diag(ones(nstate,1)) zeros(nstate)];
  
%   size(RR)
%   size(zeros(nstate,nshock))

RR = [RR; zeros(nstate,nshock)];

% compute likelihood with Kalman filter by Durbin & Koopman (2012, p85) 

t = 1;
while t <= nobs
   
   At1 = At;
   Pt1 = Pt;
   
   % Forecasting
   alphahat = TT*At1;
   
%    size(TT)
%    
%    size(Pt1)
%    
%    size(RR*QQ*RR')
   
   Phat = TT*Pt1*TT' + RR*QQ*RR';
   yhat = ZZ*alphahat + DD;
   nu   = yy(t,:) - yhat';
   
   Ft   = ZZ*Phat*ZZ' + HH + ZZ*RR*VV + (ZZ*RR*VV)';
   Ft   = 0.5*(Ft + Ft');
   
   K_g = TT*Phat*ZZ'*inv(Ft);   % Kalman Gain
  
   at_save(t,:) = alphahat;
   pt_save(:,(t-1)*2*nstate+1:t*2*nstate) = Phat;
   
   loglh = loglh -0.5*size(yy, 2)*log(2*pi)-0.5*log(det(Ft)) ...
           - 0.5*nu*inv(Ft)*nu';
   
   % Updating
   At = alphahat + (Phat*ZZ' + RR*VV)*inv(Ft)*nu';
   Pt = Phat - (Phat*ZZ'+RR*VV)*inv(Ft)*(Phat*ZZ'+RR*VV)';
   
   %  store
%    obsmean(t,:) =yhat'  ; %alphahat(nstate/2,:)';
   
    obsmean(t,:) = ( ZZ*alphahat+DD)';
   
   obsvar(t,:)  = diag(Ft)';
   nu_save(:,t) = nu';
   Ft_save(:,(t-1)*nseries+1:t*nseries) = Ft;
   Kg_save(:,(t-1)*nseries+1:t*nseries) = K_g;
   
   t = t+1;
end  

% log likelihood
rloglh = real(loglh);

% switch smooth_type
%     case 'state_smooth'
% state smoothing by Durbin & Koopman (2001, p73-75)
%                 by Durbin & Koopman (2012, p91)

r_t = zeros(nstate*2,1);
N_t = zeros(2*nstate,2*nstate); 

 for t = nobs:-1:1 

  nu = nu_save(:,t);
  Ft =  Ft_save(:,(t-1)*nseries+1:t*nseries);
  K_g = Kg_save(:,(t-1)*nseries+1:t*nseries);
  Phat = pt_save(:,(t-1)*2*nstate+1:t*2*nstate);
  L_t = TT- K_g*ZZ; 
  
  r_t = ZZ'*inv(Ft)*nu + L_t'*r_t;     % Eq.(4.32)
  N_t = ZZ'*inv(Ft)*ZZ + L_t'*N_t*L_t;
  alpha_t = at_save(t,:)' + Phat * r_t;     
  V_t   = Phat - Phat*N_t*Phat; 
  
   smooth1(t,:) =  (ZZ*alpha_t + DD)';  
      
   leverage(t,:)= (ZZ_leverage*alpha_t )';
  
 end


%     case 'disturbance_smooth' 
% disturbance smoothing by Durbin & Koopman (2001, p??)
%                       by Durbin & Koopman (2012, p96)

r_t = zeros(nstate*2,1);
N_t = zeros(2*nstate,2*nstate); 
eta = zeros(nshock,nobs);

for t = nobs:-1:1 

  nu = nu_save(:,t);
  Ft =  Ft_save(:,(t-1)*nseries+1:t*nseries);
  K_g = Kg_save(:,(t-1)*nseries+1:t*nseries);
  Phat = pt_save(:,(t-1)*2*nstate+1:t*2*nstate);
  L_t = TT- K_g*ZZ; 
  eta(:,t) = QQ*RR'*r_t;
  
  r_t = ZZ'*inv(Ft)*nu + L_t'*r_t;     % Eq.(4.32)
  N_t = ZZ'*inv(Ft)*ZZ + L_t'*N_t*L_t; 
  
end

st1 = a1;  %zeros(80,1);
 
 alpha_t = st1 + p1*r_t;  % initialize
%   alpha_t = st1;

for t = 1:nobs
   smooth2(t,:) = (ZZ * alpha_t(1:80))' + DD';
  alpha_t = TT * alpha_t + RR * eta(:,t);  % Eq.(4.85)
end;


% /*   historical decomposition
% */

nvar    = nseries;

for  index = 1:nshock;
%    st1 = zeros(40,1);
   r1 = zeros(40,1);
   r1(index) = r_t(index);  
   
%    st1 = a1 + p1*r_t;  %QQ(index,index)/sum(diag(QQ)) ;  % initialize    
    st1 = a1;                                             % initialize   

 for i =1:nobs-1;   
      st_save(i,1+(index-1)*nvar:index*nvar ) = (ZZ*st1)';
    
     shock_tmp = zeros(nshock,1);
     shock_tmp(index) = eta(index,i);     
     st = TT*st1 + RR*shock_tmp; 
     st1 = st;
 end;
 
end;

% end %( end switch)

%-----------------------------
function [omega] = createcov(para)

npara = max(size(para));
omega = zeros(npara, npara);
for i = 1:npara
  omega(i, i) = para(i)^2;
end
