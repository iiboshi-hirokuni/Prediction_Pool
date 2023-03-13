
% cd X:\jk03349\GK_JME_v2

clear all

dum1 = zeros(33,1);
dum2 = zeros(33,1);

dum1([29]) = 1;
dum2([29]) = 1;

%In_ldÇ™å¥àˆÅHÅH

%dum(2) = 1;
%dum(5) = 1;
%dum(16) = 1;

%dum(18) = 1;
%dum(25) = 1;
%dum(29) = 1;

syms betta sig hh varphi zetta alfa G_over_Y eta_i epsilon ...
     gam gam_P kappa_pi kappa_y rho_i rho_ksi sigma_ksi ...
     rho_a sigma_a rho_g sigma_g sigma_i kappa tau chi b delta_c G_ss I_ss

syms Y Ym K Keff L I C G Q varrho Lambda Rk R Pm U D X F ...
     Z i w VMPK delta In Welf a ksi g infl inflstar

syms Y_lg Ym_lg K_lg Keff_lg L_lg I_lg C_lg G_lg Q_lg ...
     varrho_lg Lambda_lg Rk_lg R_lg Pm_lg U_lg D_lg X_lg F_lg ...
     Z_lg i_lg w_lg VMPK_lg delta_lg In_lg Welf_lg a_lg ksi_lg ...
     g_lg infl_lg inflstar_lg 

syms C_ld Lambda_ld In_ld Rk_ld infl_ld F_ld Z_ld ...
     eta_C eta_Lambda eta_In eta_Rk eta_R eta_infl eta_F eta_Z ...
     C_lda Lambda_lda In_lda Rk_lda infl_lda F_lda Z_lda 

syms e_a e_ksi e_g e_i

X_ss = epsilon/(epsilon-1);

% Endogenous Variables
Y_t = [Y; Ym; K; L; I; C; G; Q; varrho; Lambda; Rk; R; Pm; U; D; X; F; ...
       Z; i; delta; In; a; ksi; g; infl; inflstar; ...
       C_ld; Lambda_ld; Rk_ld; infl_ld; F_ld; Z_ld];

% Endogenous Variables(lagged)
Y_ta = [Y_lg; Ym_lg; K_lg; L_lg; I_lg; C_lg; G_lg; Q_lg;...
        varrho_lg; Lambda_lg; Rk_lg; R_lg; Pm_lg; U_lg; D_lg; X_lg; F_lg; ...
        Z_lg; i_lg; delta_lg; In_lg; a_lg; ksi_lg; ...
        g_lg; infl_lg; inflstar_lg; ...
        C_lda; Lambda_lda; Rk_lda; infl_lda; F_lda; Z_lda];

% Forecast Errors
Eta_t = [eta_C; eta_Lambda; eta_Rk; eta_R; eta_infl; eta_F; eta_Z];

% Exogenous Shock Variables
Epsilon_t = [e_a; e_ksi; e_g; e_i];

%Household
%1. Marginal utility of consumption
Eq1 = -varrho + (C-hh*C_lg)^(-sig)-betta*hh*(C_ld-hh*C)^(-sig);

%2. Euler equation
Eq2 = -betta*R*Lambda_ld + 1.0;

%3. Stochastic discount rate
Eq3 = - Lambda + varrho/varrho_lg;

%4. Arbitrage
Eq4 = -betta*Lambda_ld*Rk_ld + betta*Lambda_ld*R;

%5. Labor market equilibrium
Eq5 = -chi*L^varphi + varrho*Pm*(1-alfa)*Ym/L;

%Final goods producer
%6. Return to capital
Eq6 = -Rk + (Pm*alfa*Ym/K_lg+ksi*(Q-delta))/Q_lg;

%7. Production function
Eq7 = -Ym + a*(ksi*U*K_lg)^alfa*L^(1-alfa);

%Capital Goods Producer
%8. Optimal investment decision
Eq8 = -Q + 1+eta_i/2*((log(In)+I_ss)/(log(In_lg)+I_ss)-1)^2+...
    eta_i*((log(In)+I_ss)/(log(In_lg)+I_ss)-1)*(log(In)+I_ss)/(log(In_lg)+I_ss)-...
    betta*Lambda_ld*eta_i*((log(1.0)+I_ss)/(log(In)+I_ss)-1)...
    *((log(1.0)+I_ss)/(log(In)+I_ss))^2;

%9. Depreciation rate
Eq9 = -delta + delta_c+b/(1+zetta)*U^(1+zetta);

%10. Optimal capacity utilization rate
Eq10 = -Pm*alfa*Ym/U + b*U^zetta*ksi*K_lg;

%11. Net investment
Eq11 = -log(In) + I-delta*ksi*K_lg;

%12. Capital accumulation equation
Eq12 = -K + ksi*K_lg+log(In); 

%13. Government consumption
Eq13 = -G + G_ss*g;

%Equilibrium
%14. Aggregate resource constraint
Eq14 = -Y + C+G+I+eta_i/2*((log(In)+I_ss)/(log(In_lg)+I_ss)-1)^2*(log(In)+I_ss);

%15. Wholesale, retail output
Eq15 = -Ym + Y*D;

%16. Price dispersion
Eq16 = -D + gam*D_lg*infl_lg^(-gam_P*epsilon)*infl^epsilon+...
      (1-gam)*((1-gam*infl_lg^(gam_P*(1-gam))*...
      infl^(gam-1))/(1-gam))^(-epsilon/(1-gam));

%17. Markup
Eq17 = -X + 1/Pm;

%18. Optimal price choice
Eq18 = -F + Y*Pm+betta*gam*Lambda_ld*infl_ld^epsilon*...
       infl^(-epsilon*gam_P)*F_ld;

%19.
Eq19 = -Z + Y+betta*gam*Lambda_ld*infl_ld^(epsilon-1)*...
       infl^(gam_P*(1-epsilon))*Z_ld;

%20. Optimal price choice
Eq20 = -inflstar + epsilon/(epsilon-1)*F/Z*infl;

%21. Price index
Eq21 = -infl^(1-epsilon) + gam*infl_lg^(gam_P*(1-epsilon))+...
      (1-gam)*(inflstar)^(1-epsilon);

%22. Fisher equation
%Eq22 = -log(i) + log(R) + log(infl_ld);
Eq22 = -i + R*infl_ld;

%23. Interest rate rule
%Eq23 = -log(i) + rho_i*log(i_lg) + (1-rho_i)*(log(1/betta) + kappa_pi*log(infl) + ...
%    kappa_y*log(X/X_ss)) + e_i;
%Eq23 = -log(i) + rho_i*log(i_lg) + (1-rho_i)*(kappa_pi*log(infl) + ...
%       kappa_y*log(X)) + e_i;

Eq23 = -log(i) + rho_i*log(i_lg) + (1-rho_i)*(log(1/betta) + kappa_pi*log(infl) + ...
    kappa_y*log(X/X_ss)) + e_i;
    
%Shocks
%24. TFP shock
Eq24 = -log(a) + rho_a*log(a_lg)-e_a;

%25. Capital quality shock
Eq25 = -log(ksi) + rho_ksi*log(ksi_lg)-e_ksi;

%26. Government consumption shock
Eq26 = -log(g) + rho_g*log(g_lg)-e_g;

EqF1 = -C + C_lda + eta_C;
EqF2 = -Lambda + Lambda_lda + eta_Lambda;
%EqF3 = -In + In_lda + eta_In;
%EqF4 = -R + R_lda + eta_R;
EqF5 = -Rk + Rk_lda + eta_Rk;
EqF6 = -infl + infl_lda + eta_infl;
EqF7 = -F + F_lda+ eta_F;
EqF8 = -Z + Z_lda+ eta_Z;

System_of_Eq = [Eq1; Eq2; Eq3; Eq4; Eq5;...
             Eq6; Eq7; Eq8; Eq9; Eq10;...
             Eq11; Eq12; Eq13; Eq14; Eq15;...
             Eq16; Eq17; Eq18; Eq19; Eq20;...
             Eq21; Eq22; Eq23; Eq24; Eq25; Eq26;...
             EqF1; EqF2; EqF5; EqF6; EqF7 ;EqF8];



neq  = length(System_of_Eq);
CC = zeros(neq,1);
GAM0j = -jacobian(System_of_Eq, Y_t);
GAM1j = jacobian(System_of_Eq, Y_ta);
PSI0j = jacobian(System_of_Eq, Epsilon_t);
PPIj  = jacobian(System_of_Eq, Eta_t);

Z0 = [Y; K; L; I; C; R; i; infl; inflstar; a; ksi; g];
ZZ = jacobian(Z0, Y_t); 

%---------------------------
% Solve the model by gensys
%---------------------------

addpath('./gensys');

% find steady state value
[Yss Kss Lss Css Iss varrhoss Pmtss Fss Zss] = findSS();

betta=0.99;
sig=1.0;
hh=0.815;
varphi=0.276;
zetta=7.2;
alfa=0.33;
G_over_Y=0.2;
eta_i=1.728;
epsilon=4.167;
gam=0.779;
gam_P=0.241;
kappa_pi=1.5;
kappa_y=-0.125;
rho_i=0.0;
rho_ksi=0.66;
sigma_ksi=0.05;
rho_a=0.95;
sigma_a=0.01;
rho_g=0.95;
sigma_g=0.01;
sigma_i=0.01;
kappa=10.0;
tau=0.001;
chi=3.47600242;
b=0.03510101;
delta_c=0.02071939;
G_Y = 0.2;

I_ss = Iss;
G_ss = Yss*G_Y;

% assign the ss value to each variables
Y = Yss;
Ym = Yss;
K = Kss;
L = Lss;
I = Iss;
C = Css;
G = G_ss;
Q = 1.0;
varrho = varrhoss;
Lambda = 1.0;
Rk = 1.0/betta;
R = 1.0/betta;
Pm = Pmtss;
U = 1.0;
D = 1.0;
X = 1.0/Pmtss;
F = Fss;
Z = Zss;
i = 1.0/betta;
delta = 0.025;
In = 1.0;
a = 1.0;
ksi = 1.0;
g = 1.0;
infl = 1.0;
inflstar = 1.0;
C_ld = C;
Lambda_ld = Lambda;
In_ld = In;
Rk_ld = Rk;
infl_ld = infl;
F_ld = F;
Z_ld = Z;

Y_lg = Y;
Ym_lg = Ym;
K_lg = K;
L_lg = L;
I_lg = I;
C_lg = C;
G_lg = G;
Q_lg = Q;
varrho_lg = varrho;
Lambda_lg = Lambda;
Rk_lg = Rk;
R_lg = R;
Pm_lg = Pm;
U_lg = U;
D_lg = D;
X_lg = X;
F_lg = F;
Z_lg = Z;
i_lg = i;
delta_lg = delta;
In_lg = In;
a_lg = a;
ksi_lg = ksi;
g_lg = g;
infl_lg = infl;
inflstar_lg = inflstar;
C_lda = C;
Lambda_lda = Lambda;
In_lda = In;
Rk_lda = Rk;
infl_lda = infl;
F_lda = F;
Z_lda = Z;

Ystar0 = eval(Y_t);
Ystar = zeros(neq, neq);

for j = 1:neq
  Ystar(j,:) = Ystar0;
end

GAM0 = eval(GAM0j).*Ystar;
GAM1 = eval(GAM1j).*Ystar;
PSI0 = eval(PSI0j);
PPI  = eval(PPIj);

[T1,TC,T0,fmat,fwt,ywt,gev,RC,loose] = gensys(GAM0,GAM1,CC,PSI0,PPI,1);

RC

%-------------------
% Impulse Response
%-------------------
%
%%Y K L I C R i infl inflstar
%
%
%nvar = size(Z0, 1);
%nirf = 40; 
%
%
%nshock = size(Epsilon_t, 1);
%sig_chol = eye(nshock); 
%
%titlestr = {'e_a :'; 'e_ksi :'; 'e_g :'; 'e_i :'};
%ystr = {'Y'; 'K'; 'L'; 'I'; 'C'; 'R'; 'i'; 'infl'; 'inflstar'; ...
%        'a :'; 'ksi :'; 'g :'};
%
%for sh_ind = 1:nshock
%    impact = sig_chol(:,sh_ind);
%    yyirf  = zeros(nirf,nvar);
%    s = T0*impact;
%    yyirf(1,:) = ZZ*s;
%    
%    for t = 2:nirf
%        ss = T1*s;     
%        yyirf(t,:) = (ZZ*ss)';
%        s = ss;
%    end
%    for j = 1:12
%      subplot(4, 3, j)
%      plot(1:nirf, yyirf(:,j),'b')
%      title(strcat(titlestr(sh_ind),ystr(j)))
%    end
%  w = waitforbuttonpress;
%end

%Some extra variables for convenience
%26. Effective capital
%exp(Keff)   =   exp(ksi)*exp(K(-1));
%
%%27. Wages
%exp(w)      =   exp(Pm)*(1-alfa)*exp(Y)/exp(L);
%
%%28. Marginal value product of capital
%exp(VMPK)   =   exp(Pm)*alfa*exp(Y)/(exp(ksi)*exp(K(-1)));
%
%%29. Welfare
%Welf   =   log(exp(C)-hh*exp(C(-1)))-chi*exp(L)^(1+varphi)/(1+varphi)+betta*Welf(+1);

%Y_t([2,5,16,18,25,29])
%Ym I X Z infl In_ld



% 26
%Y_t([13,27,31,32,33])
%Pm C_ld infl_ld F_ld Z_ld

%Y_t([27,31,32,33])
%C_ld infl_ld F_ld Z_ld
%
%dum1(31) = 1;
%dum2(13) = 1;


%dum1([6,17,18,26,31,32,33]) = 1;
%dum2([13,18,19,20,31,32,33]) = 1;

%dum1([6,17,18,26,31,32,33,31,22,23,27]) = 1;
%dum2([13,18,19,20,31,32,33,31,24,25,27]) = 1;

%dum1([12,8,1,3,10,9,11,4,5,14,20,21]) = 1; %2,7 13
%dum2([2,8,1,4,3,5,6,7,9,10,11,12]) = 1;

%dum(8) = 1;
%dum(32) = 1;
%dum(33) = 1;