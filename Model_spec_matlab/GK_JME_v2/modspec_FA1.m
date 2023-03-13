%=============================
% GK FA MODEL (JME 2011)
%=============================

% cd X:\jk03349\GK_JME_v2

clear all

syms betta sig hh varphi zetta alfa G_over_Y eta_i epsilon ...
     gam gam_P kappa_pi kappa_y rho_i rho_ksi sigma_ksi ...
     rho_a sigma_a rho_g sigma_g sigma_i kappa tau chi b delta_c G_ss I_ss ...
     lambda omega theta

syms Y Ym K L I C G Q varrho Lambda Rk R Pm U D X F ...
     Z i delta In a ksi g infl inflstar ...
     nu eta phi z x N Ne Nn prem

syms Y_lg Ym_lg K_lg L_lg I_lg C_lg G_lg Q_lg varrho_lg ...
     Lambda_lg Rk_lg R_lg Pm_lg U_lg D_lg X_lg F_lg ...
     Z_lg i_lg delta_lg In_lg a_lg ksi_lg g_lg infl_lg inflstar_lg ...
     nu_lg eta_lg phi_lg z_lg x_lg N_lg Ne_lg Nn_lg prem_lg

syms C_ld Lambda_ld In_ld Rk_ld infl_ld F_ld Z_ld ...
     nu_ld eta_ld z_ld x_ld ...
     C_lda Lambda_lda In_lda Rk_lda infl_lda F_lda Z_lda ...
     nu_lda eta_lda z_lda x_lda ...
     eta_C eta_Lambda eta_In eta_Rk eta_infl eta_F eta_Z ...
     eta_nu eta_eta eta_z eta_x ...
     e_a e_ksi e_g e_i e_Ne
     
% Endogenous Variables
Y_t = [Y; Ym; K; L; I; C; G; Q; varrho; Lambda; Rk; R; Pm; U; D; X; F; ...
       Z; i; delta; In; a; ksi; g; infl; inflstar; ...
       nu; eta; phi; z; x; N; Ne; Nn; prem; ...
       C_ld; Lambda_ld; In_ld; Rk_ld; infl_ld; F_ld; Z_ld; ...
       nu_ld; eta_ld; z_ld; x_ld];

% Endogenous Variables(lagged)
Y_ta = [Y_lg; Ym_lg; K_lg; L_lg; I_lg; C_lg; G_lg; Q_lg; varrho_lg; ...
       Lambda_lg; Rk_lg; R_lg; Pm_lg; U_lg; D_lg; X_lg; F_lg; ...
       Z_lg; i_lg; delta_lg; In_lg; a_lg; ksi_lg; g_lg; infl_lg; inflstar_lg; ...
       nu_lg; eta_lg; phi_lg; z_lg; x_lg; N_lg; Ne_lg; Nn_lg; prem_lg; ...
       C_lda; Lambda_lda; In_lda; Rk_lda; infl_lda; F_lda; Z_lda; ...
       nu_lda; eta_lda; z_lda; x_lda];
       
% Forecast Errors
Eta_t = [eta_C; eta_Lambda; eta_In; eta_Rk; eta_infl; eta_F; eta_Z; ...
         eta_nu; eta_eta; eta_z; eta_x];

% Exogenous Shock Variables
Epsilon_t = [e_a; e_ksi; e_g; e_i; e_Ne];

% Home household
% 1. Marginal utility of consumption
Eq1 = -exp(varrho) + (exp(C)-hh*exp(C_lg))^(-sig)-betta*hh*(exp(C_ld)-hh*exp(C))^(-sig);

% 2. Euler equation
Eq2 = -betta*exp(R)*exp(Lambda_ld) + 1;

% 3. Stochastic discount rate
Eq3 = -exp(Lambda) + exp(varrho)/exp(varrho_lg);

% 4. Labor market equilibrium
Eq4 = -chi*exp(L)^varphi + exp(varrho)*exp(Pm)*(1-alfa)*exp(Y)/exp(L);

% Financial Intermediaries
% 5. Value of banks' capital
Eq5 = -(exp(nu)-1) + (1-theta)*betta*exp(Lambda_ld)*(exp(Rk_ld)-exp(R))+betta*exp(Lambda_ld)*theta*exp(x_ld)*(exp(nu_ld)-1);

% 6. Value of banks' net wealth
Eq6 = -exp(eta) + (1-theta)+betta*exp(Lambda(+1))*theta*exp(z(+1))*exp(eta(+1));

% 7. Optimal leverage
Eq7 = -exp(phi) + exp(eta)/(lambda-(exp(nu)-1));

% 8. Growth rate of banks' capital
Eq8 = -exp(z) + (exp(Rk)-exp(R_lg))*exp(phi_lg)+exp(R_lg);

% 9. Growth rate of banks' net wealth
Eq9 = -exp(x) + (exp(phi)/exp(phi_lg))*exp(z);

% Aggregate capital, net worth
% 10. Aggregate capital
Eq10 = -exp(Q)*exp(K) + exp(phi)*exp(N);

% 11. Banks' net worth
Eq11 = -exp(N) + exp(Ne)+exp(Nn);

% 12. Existing banks' net worth accumulation
Eq12 = -exp(Ne) + theta*exp(z)*exp(N_lg)*exp(-e_Ne);

% 13. New banks' net worth
Eq13 = -exp(Nn) + omega*exp(Q)*exp(ksi)*exp(K_lg);

% Final goods producer
% 14. Return to capital
Eq14 = -exp(Rk) + (exp(Pm)*alfa*exp(Ym)/exp(K_lg)+exp(ksi)*(exp(Q)-exp(delta)))/exp(Q_lg);

% 15. Production function
Eq15 = -exp(Ym) + exp(a)*(exp(ksi)*exp(U)*exp(K_lg))^alfa*exp(L)^(1-alfa);

% Capital Goods Producer
% 16. Optimal investment decision
Eq16 = -exp(Q) + 1+eta_i/2*((In+I_ss)/(In_lg+I_ss)-1)^2+eta_i*((In+I_ss)/(In_lg+I_ss)-1)*(In+I_ss)/(In_lg+I_ss)-betta*exp(Lambda_ld)*eta_i*((In_ld+I_ss)/(In+I_ss)-1)*((In_ld+I_ss)/(In+I_ss))^2;

% 17. Depreciation rate
Eq17 = -exp(delta) + delta_c+b/(1+zetta)*exp(U)^(1+zetta);

% 18. Optimal capacity utilization rate
Eq18 = -exp(Pm)*alfa*exp(Ym)/exp(U) + b*exp(U)^zetta*exp(ksi)*exp(K_lg);

% 19. Net investment
Eq19 = -In + exp(I)-exp(delta)*exp(ksi)*exp(K_lg);

% 20. Capital accumulation equation
Eq20 = -exp(K) + exp(ksi)*exp(K_lg)+In; 

% 21. Government consumption
Eq21 = -exp(G) + G_ss*exp(g);

% Equilibrium
% 22. Aggregate resource constraint
Eq22 = -exp(Y) + exp(C)+exp(G)+exp(I)+eta_i/2*((In+I_ss)/(In_lg+I_ss)-1)^2*(In+I_ss);

% 23. Wholesale, retail output
Eq23 = -exp(Ym) + exp(Y)*exp(D);

% 24. Price dispersion
Eq24 = -exp(D) + gam*exp(D_lg)*exp(infl_lg)^(-gam_P*epsilon)*exp(infl)^epsilon+(1-gam)*((1-gam*exp(infl_lg)^(gam_P*(1-gam))*exp(infl)^(gam-1))/(1-gam))^(-epsilon/(1-gam));

% 25. Markup
Eq25 = -exp(X) + 1/exp(Pm);

% 26. Optimal price choice
Eq26 = -exp(F) + exp(Y)*exp(Pm)+betta*gam*exp(Lambda_ld)*exp(infl_ld)^epsilon*(exp(infl))^(-epsilon*gam_P)*exp(F_ld);

% 27.
Eq27 = -exp(Z) + exp(Y)+betta*gam*exp(Lambda_ld)*exp(infl_ld)^(epsilon-1)*exp(infl)^(gam_P*(1-epsilon))*exp(Z_ld);

% 28. Optimal price choice
Eq28 = -exp(inflstar) + epsilon/(epsilon-1)*exp(F)/exp(Z)*exp(infl);

% 29. Price index
Eq29 = -(exp(infl))^(1-epsilon) + gam*exp(infl_lg)^(gam_P*(1-epsilon))+(1-gam)*(exp(inflstar))^(1-epsilon);

% 30. Fisher equation
Eq30 = -exp(i) + exp(R)*exp(infl_ld);

% 31. Interest rate rule
Eq31 = -exp(i) + exp(i_lg)^rho_i*((1/betta)*exp(infl)^kappa_pi*(exp(X)/(epsilon/(epsilon-1)))^(kappa_y))^(1-rho_i)*exp(e_i);

% Shocks
% 32. TFP shock
Eq32 = -a + rho_a*a_lg-e_a;

% 33. Capital quality shock
Eq33 = -ksi + rho_ksi*ksi_lg-e_ksi;

% 34. Government consumption shock
Eq34 = -g + rho_g*g_lg-e_g;

% 39. Premium
Eq35 = -exp(prem) + exp(Rk_ld)/exp(R);

EqF1 = -exp(C) + exp(C_lda)*exp(eta_C);
EqF2 = -exp(Lambda) + exp(Lambda_lda)*exp(eta_Lambda);
EqF3 = -exp(In) + exp(In_lda)*exp(eta_In);
EqF4 = -exp(Rk) + exp(Rk_lda)*exp(eta_Rk);
EqF5 = -exp(infl) + exp(infl_lda)*exp(eta_infl);
EqF6 = -exp(F) + exp(F_lda)*exp(eta_F);
EqF7 = -exp(Z) + exp(Z_lda)*exp(eta_Z);

EqF8 = -exp(nu) + exp(nu_lda)*exp(eta_nu);
EqF9 = -exp(eta) + exp(eta_lda)*exp(eta_eta);
EqF10 = -exp(z) + exp(z_lda)*exp(eta_z);
EqF11 = -exp(x) + exp(x_lda)*exp(eta_x);

System_of_Eq = [Eq1; Eq2; Eq3; Eq4; Eq5;...
             Eq6; Eq7; Eq8; Eq9; Eq10;...
             Eq11; Eq12; Eq13; Eq14; Eq15;...
             Eq16; Eq17; Eq18; Eq19; Eq20;...
             Eq21; Eq22; Eq23; Eq24; Eq25;...
             Eq26; Eq27; Eq28; Eq29; Eq30;...
             Eq31; Eq32; Eq33; Eq34; Eq35;...
             EqF1; EqF2; EqF3; EqF4; EqF5;...
             EqF6; EqF7; EqF8; EqF9; EqF10; EqF11];

neq  = length(System_of_Eq);
CC = zeros(neq,1);
GAM0j = -jacobian(System_of_Eq, Y_t);
GAM1j = jacobian(System_of_Eq, Y_ta);
PSI0j = jacobian(System_of_Eq, Epsilon_t);
PPIj  = jacobian(System_of_Eq, Eta_t);

Z0 = [ksi; R; prem; Y; C; I; K; L; Q; N; infl; i];
ZZ = jacobian(Z0, Y_t); 

%Some extra variables for convenience
%35. Effective capital
%exp(Keff)   =   exp(ksi)*exp(K(-1));
%
%36. Wages
%exp(w)      =   exp(Pm)*(1-alfa)*exp(Y)/exp(L);
%
%37. Marginal value product of capital
%exp(VMPK)   =   exp(Pm)*alfa*exp(Y)/(exp(ksi)*exp(K(-1)));
%
%38. Welfare
%Welf   =   log(exp(C)-hh*exp(C(-1)))-chi*exp(L)^(1+varphi)/(1+varphi)+betta*Welf(+1);

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

%lambda = 0.381;
omega = 0.002;
theta = 0.972;

rho_i=0.0;
rho_ksi=0.66;
rho_a=0.95;
rho_g=0.95;

chi=3.47600242;
b=0.03510101;
delta_c=0.02071939;

I_ss = Iss;
G_ss = Yss*G_over_Y;

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

Nn_ss = omega*Kss;
Nn = log(Nn_ss);
Ne_ss = omega*theta*Kss/(betta-theta);
Ne = log(Ne_ss);
N = log(Ne_ss+Nn_ss);
lambda = (Ne_ss+Nn_ss)/Kss;

nu =  log(1.0);
eta = log(1.0);
phi = log(1.0/lambda);
z = log(1.0/betta);
x = log(1.0/betta);



prem = log(1.0);

C_ld = C;
Lambda_ld = Lambda;
In_ld = In;
Rk_ld = Rk;
infl_ld = infl;
F_ld = F;
Z_ld = Z;
nu_ld = nu;
eta_ld = eta;
z_ld = z;
x_ld = x;

e_a = log(1.0);
e_ksi = log(1.0);
e_g = log(1.0);
e_i = log(1.0);
e_Ne = log(1.0);

eta_C = log(1.0);
eta_Lambda = log(1.0);
eta_In = log(1.0);
eta_Rk = log(1.0);
eta_infl = log(1.0);
eta_F = log(1.0);
eta_Z = log(1.0);
eta_nu = log(1.0);
eta_eta = log(1.0);
eta_z = log(1.0);
eta_x = log(1.0);

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

nu_lg = nu;
eta_lg = eta;
phi_lg = phi;
z_lg = z;
x_lg = x;
Nn_lg = Nn;
N_lg = N;
Ne_lg = Ne;
prem_lg = prem;

C_lda = C;
Lambda_lda = Lambda;
In_lda = In;
Rk_lda = Rk;
infl_lda = infl;
F_lda = F;
Z_lda = Z;
nu_lda = nu;
eta_lda = eta;
z_lda = z;
x_lda = x;

% chek if all f evaluated at ss value are zero
check = eval(System_of_Eq)

GAM0 = eval(GAM0j);
GAM1 = eval(GAM1j);
PSI0 = eval(PSI0j);
PPI  = eval(PPIj);

[T1,TC,T0,fmat,fwt,ywt,gev,RC,loose] = gensys(GAM0,GAM1,CC,PSI0,PPI,1);
RC

%-------------------
% Impulse Response
%-------------------

nirf = 40; 

nshock = size(Epsilon_t, 1);
sig_chol = diag([0.01,0.05,0.01,0.01,0.01])*100; 

titlestr = {'eps a :'; 'eps xi :'; 'eps g :'; 'eps i :'};

ystr = {'xi'; 'R'; 'prem';'Y'; 'C'; 'I'; 'K'; 'L'; 'Q'; 'N'; 'infl'; 'i'};
index1 = [2,1];
nvar = size(Z0, 1);

for jj = 1:length(index1)
    
     figure(jj)
    
    sh_ind = index1(jj);
    impact = sig_chol(:,sh_ind);
    yyirf  = zeros(nirf,nvar);
    s = T0*impact;
    yyirf(1,:) = ZZ*s;
    for t = 2:nirf
        ss = T1*s;     
        yyirf(t,:) = (ZZ*ss)';
        s = ss;
    end
    for j = 1:12
      subplot(4, 3, j)
      plot(1:nirf, yyirf(:,j),'b')
      title(strcat(titlestr(sh_ind),ystr(j)))
    end
%     w = waitforbuttonpress;
end

