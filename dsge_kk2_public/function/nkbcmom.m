function [GAM0_c, GAM0_i, GAM0_q, GAM0_l, GAM0_w, GAM0_z, GAM0_p, GAM0_g, GAM0_m,...
GAM1_c, GAM1_i, GAM1_q, GAM1_l, GAM1_w, GAM1_z, GAM1_p, GAM1_g, GAM1_m] = nkbcmom(para, nn)

% Compute Covariance Matrix 
% and its Chol decomposition

sig_chol = eye(9);
sig_chol(1,1) = para(19);   % consumption shock
sig_chol(2,2) = para(20);   % investment  shock
sig_chol(3,3) = para(21);   % equity premium shock
sig_chol(4,4) = para(22);   % labor shock
sig_chol(5,5) = para(23);   % wage mark-up shock
sig_chol(6,6) = para(24);   % productivity shock
sig_chol(7,7) = para(25);   % price mark-up shock
sig_chol(8,8) = para(26);   % fiscal shock
sig_chol(9,9) = para(27);   % monetary policy shock

% solve the DSGE model

[TT,TC,T0,RC] = dsgesolv(para);
nstate = size(TT,2);

RR = T0;

ZZ      = zeros(7, size(TC,1));
ZZ(1,1) = 1;  % output gap
ZZ(2,2) = 4;  % inflation
ZZ(3,3) = 1;  % wage
ZZ(4,6) = 1;  % investment
ZZ(5,7) = 1;  % consumption
ZZ(6,8) = 4;  % nominal rate
ZZ(7,10) = 1; % labor

HH = zeros(7,7);
VV = zeros(9,7);

[GAM0_c, GAM1_c] = gam0gam1(sig_chol(:,1), T0, TT, ZZ, RR, VV, HH, nstate, nn);
[GAM0_i, GAM1_i] = gam0gam1(sig_chol(:,2), T0, TT, ZZ, RR, VV, HH, nstate, nn);
[GAM0_q, GAM1_q] = gam0gam1(sig_chol(:,3), T0, TT, ZZ, RR, VV, HH, nstate, nn);
[GAM0_l, GAM1_l] = gam0gam1(sig_chol(:,4), T0, TT, ZZ, RR, VV, HH, nstate, nn);
[GAM0_w, GAM1_w] = gam0gam1(sig_chol(:,5), T0, TT, ZZ, RR, VV, HH, nstate, nn);
[GAM0_z, GAM1_z] = gam0gam1(sig_chol(:,6), T0, TT, ZZ, RR, VV, HH, nstate, nn);
[GAM0_p, GAM1_p] = gam0gam1(sig_chol(:,7), T0, TT, ZZ, RR, VV, HH, nstate, nn);
[GAM0_g, GAM1_g] = gam0gam1(sig_chol(:,8), T0, TT, ZZ, RR, VV, HH, nstate, nn);
[GAM0_m, GAM1_m] = gam0gam1(sig_chol(:,9), T0, TT, ZZ, RR, VV, HH, nstate, nn);

%-----------------------------
function [GAM0ret, GAM1ret] = gam0gam1(sig_chol, T0, TT, ZZ, RR, VV, HH, nstate, nn)

QQ = sig_chol*sig_chol';
GAM0s = T0*QQ*T0' ;
AA = eye(nstate);
for t = 1:nn
   AA = TT*AA;
   G  = AA*T0*QQ*T0'*AA' ;
   GAM0s = GAM0s + G;
end

GAM0ret = ZZ*GAM0s*ZZ' + ZZ*RR*VV + (ZZ*RR*VV)' + HH;
GAM1ret = ZZ*(TT*GAM0s)*ZZ' + ZZ*(TT*RR*VV);

