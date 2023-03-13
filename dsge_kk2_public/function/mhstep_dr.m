% ===================================================
%   mhstep_dr
% ===================================================

% mhstep_dr()
% mhstep_fcn()
% simulation_smoother()
% draw_coef_g()
% draw_coef()
% draw_sig()

function mhstep_dr(sigmult,hessian,nblock,nsim,cc0,cc,paras,vst,runpath,...
    appendmode,lastGoodBlock,resupath,nchain,sysmat,casefile)

npara = length(paras);
para = zeros(npara,1);
for i = 1:npara
    para(i,1) = paras(i);
end

ndata = size(vst.data , 2);
nobs  = size(vst.data , 1);
nstate = size(sysmat.GAM0, 2);
nshock = size(sysmat.PSI0, 2);

% best to display some model identification stuff here?

sigscale = sigmult;
sigpropinv = hessian;

sigpropdim = sum(vst.pmaskinv);

[u,s,v] = svd(sigpropinv);
sigproplndet = 0.0;

for i = 1:npara
    if i <= sigpropdim
        s(i,i) = 1.0/s(i,i);
        sigproplndet = sigproplndet + log(s(i,i));
    end
end

% ----------------------------------
%   Metropolis Hastings Algorithm
%

%  Prior and Save area
%shocksum = zeros( nobs-1, nshock);
%decomp_sum = zeros( nobs-1, nshock*ndata);

% Initialize PARA NEW
Lam = make_zz(casefile, nstate);
Sig = 0.1*eye(ndata);
Psi = 0.1*ones(ndata,1);

while 1
     
    paranew = para + cc0*(sigscale*randn(npara,1));
    paranew = paranew .* vst.pmaskinv + vst.pfix .* vst.pmask;
    
    [postnew,likenew] = mhstep_fcn(paranew, vst, sysmat, Sig, Psi, Lam, casefile);
    
    propdens = -0.5*sigpropdim*log(2*pi) - 0.5*sigproplndet - ...
       0.5*sigpropdim*log(cc0*cc0) - 0.5*(paranew-para)'* ...
       sigpropinv*(paranew-para)/(cc0*cc0);
    
    %if postnew > -10000;
    %if postnew > -1E+10;
    if postnew > -1E+6;
        break;
    end
    
end
disp('chain initialized');

newnblock = nblock - lastGoodBlock;

for n = 1:newnblock
    
    parasim     = zeros(nsim,npara);
    likesim     = zeros(nsim,1);
    postsim     = zeros(nsim,1);
    rej         = zeros(nsim,1);
    propsim     = zeros(nsim,1);
    proppostsim = zeros(nsim,1);
    
    sigmasim    = zeros(nsim, ndata);
    psisim      = zeros(nsim, ndata);
    gamsim      = [];
    gam_mat_sim = zeros(nsim, ndata*nstate);
    
    state_save  = zeros(nsim, nstate*nobs);
    shock_save  = zeros(nsim, nshock*(nobs-1));
    
    if n == 1 && appendmode ~= 1
        
        newj = 2;
        likesim(1,1) = likenew;
        postsim(1,1) = postnew;
        propsim(1,1) = propdens;
        proppostsim(1,1) = postnew;
        parasim(1,:) = paranew';
        postold = postnew;
        likeold = likenew;
        paraold = paranew;
        
    elseif n == 1 && appendmode == 1
        
        newj = 1;
        postold = postnew;
        likeold = likenew;
        paraold = paranew;
        
    else
        
        newj = 1;
        
    end
    
    for j = newj:nsim
        
        % Gibs Sampling
        [statenew, shocknew] = ...
            simulation_smoother(paraold, vst.data, sysmat, Sig, Psi, Lam, casefile);
        [Lam, gam_g] = draw_Lam(vst.data, statenew, Sig, Psi, Lam, casefile);
        Sig = draw_Sig(vst.data, statenew, Psi, Lam, casefile); 
        Psi = draw_Psi(vst.data, statenew, Sig, Lam, casefile);
        
        % MH Step
        paranew = paraold + cc*(sigscale*randn(npara,1));
        paranew = paranew .* vst.pmaskinv + vst.pfix .* vst.pmask;
        
        [postold,likeold] = ...
           mhstep_fcn(paraold, vst, sysmat, Sig, Psi, Lam, casefile);
        [postnew,likenew] = ...
           mhstep_fcn(paranew, vst, sysmat, Sig, Psi, Lam, casefile);
        
        propdens = -0.5*sigpropdim*log(2*pi) - 0.5*sigproplndet - ...
            0.5*sigpropdim*log(cc*cc) - ...
            0.5*(paranew-paraold)'*sigpropinv*(paranew-paraold)/(cc*cc);
        
        propsim(j,1) = propdens;
        proppostsim(j,1) = postnew;
        
        r = min(1,exp(postnew-postold));
        
        if rand < r
            postsim(j,1) = postnew;
            likesim(j,1) = likenew;
            parasim(j,:) = paranew';
            
            paraold = paranew;
        else
            likesim(j,1) = likeold;
            postsim(j,1) = postold;
            parasim(j,:) = paraold';
            rej(j,1) = 1;
        end
        
        psisim(j,:)      = Psi';
        %gamsim(i,:)      = gam_g';  
	    gamsim           = [gamsim; gam_g'];
        gam_mat_sim(j,:) = reshape(Lam, 1, ndata*nstate);
        sigmasim(j,:)    = sqrt(diag(Sig)');
        
        state_save(j,:) = reshape(statenew, 1, nstate*nobs);
        shock_save(j,:) = reshape(shocknew', 1, nshock*(nobs-1));
        
        %if mod(j,nsim/10) == 0;
        %    disp(sprintf('block %d of %d : sim % d of %d',...
        %        (n+lastGoodBlock),nblock,j,nsim));
        %end
        
    end
    
    
    disp(sprintf('chain %d block %d of %d : sim % d of %d',...
                 nchain,(n+lastGoodBlock),nblock,j,nsim));
    disp(sprintf('Rejection rate:    %f\n',mean(rej)));
    disp(sprintf('Likelihood:        %f\n',mean(likesim)));
    disp(sprintf('Posterior:         %f\n',mean(postsim)));
    
    fname = strcat(resupath, runpath,'/-block-',num2str(n+lastGoodBlock),...
                   '_chain',num2str(nchain));
    save(fname, 'parasim', 'postsim', 'likesim', 'proppostsim', 'rej',...
         'psisim', 'gamsim', 'gam_mat_sim', 'sigmasim', 'state_save',...
         'shock_save');
end

% ======================================
% MH sampling function definition
% ======================================
function [lnpost,lnpy] = mhstep_fcn(para, v, sysmat, HH, Psi, zz, casefile)

bounds = v.trspec(:,[2,3]);

ind1 = para > bounds(:,1);
ind2 = para < bounds(:,2);

ind1 = delif(ind1,v.pmask);
ind2 = delif(ind2,v.pmask);

if all(ind1) && all(ind2)
    
    [lnpy,retcode,obserror,obsvar] = ...
        evaldsge_dr(para, v.data, sysmat, HH, Psi, zz, casefile);
    lnpy = real(lnpy);
    lnprio = priodens(para, v.pmean, v.pstdd, v.pshape);
    lnpost = lnpy + real(lnprio);
    
else
    
    lnpost = -1e6;
    lnpy = -1e20;
     
end

% ==============================================
% Simulation Smoother
% ==============================================
function [obssample, eta, st_save] = ...
    simulation_smoother(para, yy, sysmat, HH, Psi, zz1, casefile);

npara   = size(para, 1);
retcode = 0;

% solve the DSGE model

[T1,TC,TEPS,RC] = dsgesolv(para, sysmat);

nseries  = size(yy, 2);
nstate   = size(T1, 2)*2;
nshock = size(TEPS, 2);

nobs      = size(yy, 1)-1;
nvar      = size(yy, 2);
loglh     = 0;
loglhzero = -1E8;
obsmean   = zeros(nobs, nseries);
obsvar    = zeros(nobs, nseries);

nu_save   = zeros(nseries,nobs);
Ft_save   = zeros(nseries,nseries*nobs); 
Kg_save   = zeros(nstate,nseries*nobs);
At_save   = zeros(nstate,nobs);
Pt_save   = zeros(nstate,nstate*nobs);

% -------------------
% Check determinacy 
%
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
      zeros(nstate/2,nshock)];

ZZ = [zz1, (-1)*diag(Psi)*zz1];

% Check whether covariance matrix QQ is positive definite

if sum(eig(QQ) < 0) > 0
   loglh = loglhzero;
   return;
end

% We can now define the initial mean and variance for the state vector

yy_star = zeros(nobs, nvar);
for i = 1:nobs
  yy_star(i,:) = yy(i+1,:) - yy(i,:).*Psi';
end

At = make_init(yy_star, casefile, nstate);
Pt = dlyap(TT,RR*QQ*RR');

A1 = At; 
P1 = Pt;

% compute likelihood with Kalman filter

t = 1;
while t <= nobs
   
   At1 = At;
   Pt1 = Pt;
   
   % Forecasting
   alphahat = TT*At1;
   Phat = TT*Pt1*TT' + RR*QQ*RR';
   yhat = ZZ*alphahat + DD;
   nu   = yy_star(t,:) - yhat';  %
   
   Ft = ZZ*Phat*ZZ' + HH + ZZ*RR*VV + (ZZ*RR*VV)';
   Ft = 0.5*(Ft + Ft');
   
   K_g = TT*Phat*ZZ'*inv(Ft);   % Kalman Gain
   Pt_save(:,(t-1)*nstate+1:t*nstate) = Phat;
   
   %loglh = loglh -0.5*size(yy, 2)*log(2*pi)-0.5*log(det(Ft)) ...
   %        - 0.5*nu*inv(Ft)*nu';
   
   % Updating
   At = alphahat + (Phat*ZZ' + RR*VV)*inv(Ft)*nu';
   Pt = Phat - (Phat*ZZ'+RR*VV)*inv(Ft)*(Phat*ZZ'+RR*VV)';
   
   %  store
   obsmean(t,:) = yhat';
   obsvar(t,:)  = diag(Ft)';
   nu_save(:,t) = nu';
   Ft_save(:,(t-1)*nseries+1:t*nseries) = Ft;
   Kg_save(:,(t-1)*nseries+1:t*nseries) = K_g;
   
   t = t+1;
end  


% Simulation smoother  by DeJong & Shephard (1995)

r_t = zeros(nstate,1);
N_t = zeros(nstate,nstate); 
eta = zeros(nshock,nobs);
d = randn(nshock,nobs); %%%

t = nobs; 
while t > 0
    nu = nu_save(:,t);
    Ft =  Ft_save(:,(t-1)*nseries+1:t*nseries);
    K_g = Kg_save(:,(t-1)*nseries+1:t*nseries);
    L_t = TT - K_g*ZZ; 
    
    W_t = QQ * RR' * N_t * L_t;               % 4.78 
    c_t = QQ - QQ * RR' * N_t * RR * QQ;      % 4.83
    c_t = 0.5 * (c_t + c_t');   
    
    % positive definite
%    if ( det(c_t(1,1))>0 ) && (det(c_t(1:2,1:2))>0) && (det(c_t(1:3,1:3))>0)...
%        && (det(c_t(1:4,1:4))>0) && (det(c_t(1:5,1:5))>0) && (det(c_t(1:6,1:6))>0)...
%        && (det(c_t(1:7,1:7))>0) && (det(c_t(1:8,1:8))>0) && (det(c_t(1:9,1:9))>0)
%       
        e= chol(c_t)*d(:,t);
        
%    else
%        e= diag(c_t).*d(:,t);
%        disp('c_t is not positive definite'); 
%    end
    
    eta(:,t) = e + QQ*RR'*r_t;      % 4.84
    
    %  backward   
    r_t = ZZ'*inv(Ft)*nu - W_t'*inv(c_t)*e + L_t'*r_t;        %4.75
    N_t = ZZ'*inv(Ft)*ZZ + W_t'*inv(c_t)*W_t + L_t'*N_t*L_t;
    t = t-1;
end

% sampling state variables


alpha_t = A1 + P1*r_t;  % initialize
obssample = zeros(nobs+1, nstate/2);
obssample(1,:) = alpha_t(nstate/2+1:nstate);
for t = 1:nobs
   obssample(t+1,:) = (alpha_t(1:nstate/2))';
   alpha_t = TT * alpha_t + RR * eta(:,t);   % 4.85
end

% historical decomposition

st_save = zeros(nobs+1,nvar*nshock);
for index = 1:nshock
    
    st1 = (A1 + P1*r_t)*QQ(index,index)/sum(sum(diag(QQ))); % initialize
    st_save(1,1+(index-1)*nvar:index*nvar) = (zz1*st1(nstate/2+1:nstate))';
    
    for t = 1:nobs
        st_save(t+1,1+(index-1)*nvar:index*nvar) = (zz1*st1(1:nstate/2))';
        shock_tmp = zeros(nshock,1);
        shock_tmp(index) = eta(index,t);
        st = TT*st1 + RR*shock_tmp;
        st1 = st;
  end
end

% ==============================================
% createcov
% ==============================================
function omega = createcov(para)
    npara = max(size(para));
    omega = zeros(npara, npara);
    for i = 1:npara
        omega(i, i) = para(i)^2;
    end

% ==============================================
% draw_coef_g
% ==============================================
function [Lam, gam_g] = draw_Lam(yy, ss, Sig, Psi, Lam0, casefile) ;
    
    Lam   = Lam0;
	num   = size(yy, 2);    % number of observable variables
    nobs  = size(yy, 1);
    ystar = zeros(nobs-1,num);
    for i = 1:nobs-1
      ystar(i,:) = yy(i+1,:) - yy(i,:).*Psi';
    end
    
	n_index = size(ss,2);
    gam_g = [];
    
    for i = 1:num
        
        %k = casefile(i,1);
        l = casefile(i,2); 
        
        if casefile(i,4) == 2 % sensor
            
            t0 = casefile(i,10);
            R0 = casefile(i,11);
            
            xstar = ss(2:nobs,l) - Psi(i)* ss(1:nobs-1,l);
            sig2_v   = Sig(i,i);
            V        = inv(R0^(-1) + sig2_v^(-1)*xstar(:,1)'* xstar(:,1));
            gam_tmp  = V*(R0^(-1)*t0 + sig2_v^(-1)*xstar(:,1)'* ystar(:,i));
            C        = chol(V);
            gam_g0   = gam_tmp + C'*randn(1,1);
			gam_g    = [gam_g; gam_g0];
			Lam(i,l) = gam_g0;
            
    	elseif casefile(i,4) == 3 % info
		    
		    xstar = ss(2:nobs,:)' - Psi(i)* ss(1:nobs-1,:)';
			
			t0 = ones(n_index,1)*casefile(i,10);
			R0 = eye(n_index)*casefile(i,11)^(-1); 
            
			sig2_v   = Sig(i,i);           
            V        = inv(R0 + sig2_v^(-1)*xstar'* xstar);
            gam_tmp  = V*(R0*t0 + sig2_v^(-1)*xstar'* ystar(:,i)); 
			C        = chol(V);
            gam_g0   = gam_tmp + C'*randn(n_index,1);
			gam_g    = [gam_g; gam_g0];
			Lam(i,:) = gam_g0';
        end
    end


% ==============================================
% draw_coef
% ==============================================
function Psi = draw_Psi(yy, ss, Sig, Lam, casefile);
    
    num = size(yy, 2);    % number of  observable variables
    resid  = yy - ss*Lam';
    ystar = resid(2:size(resid,1),:);
    xstar = resid(1:size(resid,1)-1,:);
    tstar = size(ystar, 1);
    
    Psi0   = zeros(num,1);
    Psi = zeros(num,1);
    
    for i = 1:num
        if casefile(i,2)== 4; % nominal interest rate
            Psi(i) = 0;
        else
            t0 = casefile(i,6);
            R0 = casefile(i,7);
            
            sig2_v = Sig(i,i);
            V = inv(R0^(-1) + sig2_v^(-1)*xstar(:,i)'* xstar(:,i));
            Psi0(i) =  V*(R0^(-1)*t0 + sig2_v^(-1)*xstar(:,i)'* ystar(:,i));
            C = chol(V);
            
            accept = 0;
            while accept == 0;
                Psi(i) = Psi0(i) + C'*randn(1,1);
                if abs(Psi(i)) < 1
                    accept = 1;
                else
                    accept = 0;
                end
            end
        end
    end


% ==============================================
% draw_sig
% ==============================================
function omega = draw_Sig(yy, ss, Psi, Lam, casefile);
    
    num = size(yy, 2);    % number of  observable variables
    nobs = size(yy,1)-1;
    ystar = zeros(nobs, num);
    xstar = zeros(nobs, num);
    for i = 2:nobs+1
        ystar(i-1,:) = yy(i,:) - yy(i-1,:).*Psi';
        xstar(i-1,:) = ss(i,:)*Lam' - ss(i-1,:)*Lam'.*Psi';
    end
    resid = ystar - xstar;
    
    sig2_v = 0;
    omega = zeros(num,num);
    
    for i = 1:num
        
	    v0_ = casefile(i,8);
	    d0_ = casefile(i,9);
        
        nn = nobs + v0_;
        d = d0_ + (resid(:,i))'*(resid(:,i));
        c = rndc(nn);
        t2 = c/d;
        sig2_v = 1/t2;
        omega(i,i) = sig2_v; 
        
    end

