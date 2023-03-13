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

