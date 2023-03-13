
% cd X:\jk03349\GK_JME_v2

clear all

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
       C_ld; Lambda_ld; In_ld; Rk_ld; infl_ld; F_ld; Z_ld];

% Endogenous Variables(lagged)
Y_ta = [Y_lg; Ym_lg; K_lg; L_lg; I_lg; C_lg; G_lg; Q_lg;...
        varrho_lg; Lambda_lg; Rk_lg; R_lg; Pm_lg; U_lg; D_lg; X_lg; F_lg; ...
        Z_lg; i_lg; delta_lg; In_lg; a_lg; ksi_lg; ...
        g_lg; infl_lg; inflstar_lg; ...
        C_lda; Lambda_lda; In_lda; Rk_lda; infl_lda; F_lda; Z_lda];

% Forecast Errors
Eta_t = [eta_C; eta_Lambda; eta_Rk; eta_R; eta_infl; eta_F; eta_Z];

% Exogenous Shock Variables
Epsilon_t = [e_a; e_ksi; e_g; e_i];

%Household
%1. Marginal utility of consumption
Eq1 = -exp(varrho) + (exp(C)-hh*exp(C_lg))^(-sig)-betta*hh*(exp(C_ld)-hh*exp(C))^(-sig);

%2. Euler equation
Eq2 = -betta*exp(R)*exp(Lambda_ld) + 1;

%3. Stochastic discount rate
Eq3 = -exp(Lambda) + exp(varrho)/exp(varrho_lg);

%4. Arbitrage
Eq4 = -betta*exp(Lambda_ld)*exp(Rk_ld) + betta*exp(Lambda_ld)*exp(R);

%5. Labor market equilibrium
Eq5 = -chi*exp(L)^varphi + exp(varrho)*exp(Pm)*(1-alfa)*exp(Ym)/exp(L);

%Final goods producer
%6. Return to capital
Eq6 = -exp(Rk) + (exp(Pm)*alfa*exp(Ym)/exp(K_lg)+exp(ksi)*(exp(Q)-exp(delta)))/exp(Q_lg);

%7. Production function
Eq7 = -exp(Ym) + exp(a)*(exp(ksi)*exp(U)*exp(K_lg))^alfa*exp(L)^(1-alfa);

%Capital Goods Producer
%8. Optimal investment decision
Eq8 = -exp(Q) + 1+eta_i/2*((In+I_ss)/(In_lg+I_ss)-1)^2+eta_i*((In+I_ss)/(In_lg+I_ss)-1)*(In+I_ss)/(In_lg+I_ss)-betta*exp(Lambda_ld)*eta_i*((In_ld+I_ss)/(In+I_ss)-1)*((In_ld+I_ss)/(In+I_ss))^2;

%9. Depreciation rate
Eq9 = -exp(delta) + delta_c+b/(1+zetta)*exp(U)^(1+zetta);

%10. Optimal capacity utilization rate
Eq10 = -exp(Pm)*alfa*exp(Ym)/exp(U) + b*exp(U)^zetta*exp(ksi)*exp(K_lg);

%11. Net investment
Eq11 = -In + exp(I)-exp(delta)*exp(ksi)*exp(K_lg);

%12. Capital accumulation equation
Eq12 = -exp(K) + exp(ksi)*exp(K_lg)+In; 

%13. Government consumption
Eq13 = -exp(G) + G_ss*exp(g);

%Equilibrium
%14. Aggregate resource constraint
Eq14 = -exp(Y) + exp(C)+exp(G)+exp(I)+eta_i/2*((In+I_ss)/(In_lg+I_ss)-1)^2*(In+I_ss);

%15. Wholesale, retail output
Eq15 = -exp(Ym) + exp(Y)*exp(D);

%16. Price dispersion
Eq16 = -exp(D) + gam*exp(D_lg)*exp(infl_lg)^(-gam_P*epsilon)*exp(infl)^epsilon+(1-gam)*((1-gam*exp(infl_lg)^(gam_P*(1-gam))*exp(infl)^(gam-1))/(1-gam))^(-epsilon/(1-gam));

%17. Markup
Eq17 = -exp(X) + 1/exp(Pm);

%18. Optimal price choice
Eq18 = -exp(F) + exp(Y)*exp(Pm)+betta*gam*exp(Lambda_ld)*exp(infl_ld)^epsilon*(exp(infl))^(-epsilon*gam_P)*exp(F_ld);

%19.
Eq19 = -exp(Z) + exp(Y)+betta*gam*exp(Lambda_ld)*exp(infl_ld)^(epsilon-1)*exp(infl)^(gam_P*(1-epsilon))*exp(Z_ld);

%20. Optimal price choice
Eq20 = -exp(inflstar) + epsilon/(epsilon-1)*exp(F)/exp(Z)*exp(infl);

%21. Price index
Eq21 = -(exp(infl))^(1-epsilon) + gam*exp(infl_lg)^(gam_P*(1-epsilon))+(1-gam)*(exp(inflstar))^(1-epsilon);

%22. Fisher equation
Eq22 = -exp(i) + exp(R)*exp(infl_ld);

%23. Interest rate rule
Eq23 = -i + log(exp(i_lg)^rho_i*((1/betta)*exp(infl)^kappa_pi*(exp(X)/(epsilon/(epsilon-1)))^(kappa_y))^(1-rho_i))+e_i;

%Shocks
%23. TFP shock
Eq24 = -a + rho_a*a_lg-e_a;

%24. Capital quality shock
Eq25 = -ksi + rho_ksi*ksi_lg-e_ksi;

%25. Government consumption shock
Eq26 = -g + rho_g*g_lg-e_g;

EqF1 = -exp(C) + exp(C_lda) + eta_C;
EqF2 = -exp(Lambda) + exp(Lambda_lda) + eta_Lambda;
EqF3 = -exp(In) + exp(In_lda) + eta_In;
EqF4 = -exp(Rk) + exp(Rk_lda) + eta_Rk;
EqF5 = -exp(infl) + exp(infl_lda) + eta_infl;
EqF6 = -exp(F) + exp(F_lda)+ eta_F;
EqF7 = -exp(Z) + exp(Z_lda)+ eta_Z;

System_of_Eq = [Eq1; Eq2; Eq3; Eq4; Eq5;...
             Eq6; Eq7; Eq8; Eq9; Eq10;...
             Eq11; Eq12; Eq13; Eq14; Eq15;...
             Eq16; Eq17; Eq18; Eq19; Eq20;...
             Eq21; Eq22; Eq23; Eq24; Eq25; Eq26;...
             EqF1; EqF2; EqF3; EqF4; EqF5; EqF6 ;EqF7];

%Y_t = Y_t((dum(:,1)==0));
%Y_ta = Y_ta((dum(:,1)==0));
%System_of_Eq = System_of_Eq((dum(:,1)==0));

neq  = length(System_of_Eq);
CC = zeros(neq,1);
GAM0j = -jacobian(System_of_Eq, Y_t);
GAM1j = jacobian(System_of_Eq, Y_ta);
PSI0j = jacobian(System_of_Eq, Epsilon_t);
PPIj  = jacobian(System_of_Eq, Eta_t);

Z0 = [Y; K; L; I; C; R; i; infl; inflstar];
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
Y = log(Yss);
Ym = log(Yss);
K = log(Kss);
L = log(Lss);
I = log(Iss);
C = log(Css);
G = log(G_ss);
Q = log(1.0);
varrho = log(varrhoss);
Lambda = log(1.0);
Rk = log(1.0/betta);
R = log(1.0/betta);
Pm = log(Pmtss);
U = log(1.0);
D = log(1.0);
X = log(1.0/Pmtss);
F = log(Fss);
Z = log(Zss);
i = log(1.0/betta);
delta = log(0.025);
In = log(1.0);
a = log(1.0);
ksi = log(1.0);
g = log(1.0);
infl = log(1.0);
inflstar = log(1.0);
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
  Ystar(j,:) = exp(Ystar0);
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

%Y K L I C R i infl inflstar


nvar = size(Z0, 1);
nirf = 40; 


nshock = size(Epsilon_t, 1);
sig_chol = eye(nshock); 

titlestr = {'e_a :'; 'e_ksi :'; 'e_g :'; 'e_i :'};
ystr = {'Y'; 'K'; 'L'; 'I'; 'C'; 'R'; 'i'; 'infl'; 'inflstar'};

for sh_ind = 1:nshock
    impact = sig_chol(:,sh_ind);
    yyirf  = zeros(nirf,nvar);
    s = T0*impact;
    yyirf(1,:) = ZZ*s;
    
    for t = 2:nirf
        ss = T1*s;     
        yyirf(t,:) = (ZZ*ss)';
        s = ss;
    end
    for j = 1:9
      subplot(3, 3, j)
      plot(1:nirf, yyirf(:,j),'b')
      title(strcat(titlestr(sh_ind),ystr(j)))
    end
  w = waitforbuttonpress;
end

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

