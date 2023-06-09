%===================================================
% Modelspec.m
%   Define the Model of Kaihatsu and Kurozumi[2010]
%===================================================

para_names = {'sigma','theta','kai','inv_zeta',...
     'mu','phi_o_y','gamma_w','ksi_w', ...
     'gamma_p','ksi_p','phi_r','phi_pi','phi_y','z_star_bar','psi_bar', ...
     'eta','n_k','mu_E','r_E_bar', ...
     'rho_b','rho_g','rho_w','rho_p','rho_r','rho_nu','rho_z','rho_psi', ...
     'rho_efp','rho_nw', ...
     'sigma_b','sigma_g','sigma_w','sigma_p','sigma_r', ...
     'sigma_nu','sigma_z','sigma_psi','sigma_efp','sigma_nw'};

syms sigma theta kai inv_zeta mu phi_o_y gamma_w ksi_w ...
     gamma_p ksi_p phi_r phi_pi phi_y z_star_bar psi_bar ...
     eta n_k mu_E r_E_bar ...
     rho_b rho_g rho_w rho_p rho_r rho_nu rho_z rho_psi ...
     rho_efp rho_nw ...
     sigma_b sigma_g sigma_w sigma_p sigma_r ...
     sigma_nu sigma_z sigma_psi sigma_efp sigma_nw
     
syms z_star psi l l_bar pii pi_bar r_n r_n_bar r_E c_y i_y ...
     delta alpha lambda_w lambda_i

syms c_t lambda_t w_t b_t Et_r_E_t1 q_t n_t r_k_t u_t ...
     mc_t pi_t pi_t_lg1 pi_t_lg2 ...
     y_t i_t k_t r_n_t y_star_t ...
     z_b_t z_g_t z_w_t z_p_t z_nu_t z_r_t z_z_t z_psi_t ...
     z_efp_t z_nw_t z_star_t l_t ...
     Et_c_t1 Et_lambda_t1 Et_w_t1 Et_q_t1 Et_r_k_t1 Et_pi_t1 Et_i_t1 ...
     Et_z_b_t1 Et_z_psi_t1 Et_z_star_t1
     
syms c_ta lambda_ta w_ta b_ta Et_r_E_t1a q_ta n_ta r_k_ta u_ta ...
     mc_ta pi_ta pi_t_lg1a pi_t_lg2a ...
     y_ta i_ta k_ta r_n_ta y_star_ta ...
     z_b_ta z_g_ta z_w_ta z_p_ta z_nu_ta z_r_ta z_z_ta z_psi_ta ...
     z_efp_ta z_nw_ta z_star_ta l_ta ...
     Et_c_t1a Et_lambda_t1a Et_w_t1a Et_q_t1a Et_r_k_t1a Et_pi_t1a ...
     Et_i_t1a Et_z_b_t1a Et_z_psi_t1a Et_z_star_t1a

syms eta_c eta_lambda eta_w eta_q eta_r_k eta_pi eta_i ...
     eta_z_b eta_z_psi eta_z_star

syms epsilon_b epsilon_g epsilon_w epsilon_p epsilon_nu epsilon_r ...
     epsilon_z epsilon_psi epsilon_efp epsilon_nw
 
 syms Et_chi_t1 chi_t eta_chi Et_chi_t1a chi_ta r_E_t r_E_ta ...
      z_i_t z_i_ta eta_chi epsilon_i rho_i sigma_i

parameters = [sigma; theta; kai; inv_zeta; mu; phi_y; gamma_w; ksi_w; ...
     gamma_p; ksi_p; phi_r; phi_pi; phi_y; z_star_bar; psi_bar; ...
     eta; n_k; mu_E; r_E_bar; ...
     rho_b; rho_g; rho_w; rho_p; rho_r; rho_nu; rho_z; rho_psi; ...
     rho_efp; rho_nw; ...
     sigma_b; sigma_g; sigma_w; sigma_p; sigma_r; ...
     sigma_nu; sigma_z; sigma_psi; sigma_efp; sigma_nw];

%% Endogenous Variables (40)
Y_t = [c_t; lambda_t; w_t; b_t; r_E_t; q_t; n_t; r_k_t; u_t; ...
     mc_t; pi_t; pi_t_lg1; pi_t_lg2; ...
     y_t; i_t; k_t; r_n_t; ...
     z_b_t; z_g_t; z_w_t; z_p_t; z_nu_t; z_r_t; z_z_t; z_psi_t; ...
     z_efp_t; z_nw_t; z_star_t; l_t; ...
     Et_chi_t1; chi_t; z_i_t; ...
     Et_c_t1; Et_lambda_t1; Et_w_t1; Et_r_k_t1; Et_pi_t1; Et_i_t1; ...
     Et_z_b_t1;  Et_z_star_t1];

%% Endogenous Variables(lagged) (40)
Y_ta = [c_ta; lambda_ta; w_ta; b_ta; r_E_ta; q_ta; n_ta; r_k_ta; u_ta; ...
     mc_ta; pi_ta; pi_t_lg1a; pi_t_lg2a; ...
     y_ta; i_ta; k_ta; r_n_ta; ...
     z_b_ta; z_g_ta; z_w_ta; z_p_ta; z_nu_ta; z_r_ta; z_z_ta; z_psi_ta; ...
     z_efp_ta; z_nw_ta; z_star_ta; l_ta; ...
     Et_chi_t1a; chi_ta; z_i_ta; ...
     Et_c_t1a; Et_lambda_t1a; Et_w_t1a; Et_r_k_t1a; Et_pi_t1a; Et_i_t1a; ...
     Et_z_b_t1a;  Et_z_star_t1a];

%% Forecast Errors (9)
Eta_t = [eta_c; eta_lambda; eta_w; eta_chi; eta_r_k; eta_pi; eta_i; ...
         eta_z_b; eta_z_star];

%% Exogenous Shock Variables (11)
Epsilon_t = [epsilon_b; epsilon_g; epsilon_w; epsilon_p; epsilon_nu; ...
             epsilon_r; epsilon_z; epsilon_psi; epsilon_i; epsilon_efp; epsilon_nw];

%% Structual Equations(17)
Eq_S1 = -lambda_t - 1/(1-theta*pii/r_n)*(sigma/(1-theta/z_star)* ...
        (c_t-theta/z_star*(c_ta-z_star_t))-z_b_t) + ...
        theta*pii/r_n/(1-theta*pii/r_n)*(sigma/(1-theta/z_star)* ...
        (Et_c_t1+Et_z_star_t1-theta/z_star*c_t) - Et_z_b_t1);

Eq_S2 = -lambda_t+Et_lambda_t1-sigma*Et_z_star_t1+r_n_t-Et_pi_t1;

Eq_S3 = -w_t+w_ta-pi_t+gamma_w*pi_ta-z_star_t+z_star*pii/r_n* ...
      (Et_w_t1-w_t+Et_pi_t1-gamma_w*pi_t+Et_z_star_t1) +...
      (1-ksi_w)*(1-ksi_w*z_star*pii/r_n)/ksi_w/(1+kai*(1+lambda_w)/lambda_w)*...
      (kai*l_t-lambda_t-w_t+z_b_t) + z_w_t;

Eq_S4 = -b_t + (1+lambda_i)/(1+lambda_i-n_k)*(q_t+k_t)+...
        (1-(1+lambda_i)/(1+lambda_i-n_k))*n_t;

 Eq_S5 = -chi_t + (1-(1-delta)/r_n/psi)*r_k_t+...   
        ((1-delta)/r_n/psi)*q_t - q_ta - z_psi_t;  

Eq_S6 = -r_E_t + r_n_t - mu_E*(n_t-q_t-k_t) + z_efp_t;

Eq_S7 = -z_star/eta/r_E*n_t + (1+lambda_i)/n_k*...
        ((1-(1-delta)/r_E/psi)*r_k_t+(1-delta)/r_E/psi*q_t-q_ta-z_psi_t)-...
        ((1+lambda_i)/n_k-1)*Et_r_E_t1a+ ...
        n_ta - z_star_t + z_nw_t;

Eq_S8 = w_t + l_t - (r_k_t+u_t+k_ta-z_star_t-z_psi_t);
Eq_S9 = -u_t + mu*(r_k_t-q_t);
Eq_S10 = -mc_t + (1-alpha)*w_t + alpha*r_k_t;

Eq_S11 = -pi_t + gamma_p*pi_ta+z_star*pii/r_n*(Et_pi_t1-gamma_p*pi_t)+...
         (1-ksi_p)*(1-ksi_p*z_star*pii/r_n)/ksi_p*mc_t+z_p_t;

Eq_S12 = -y_t + (1+phi_o_y)*((1-alpha)*l_t+alpha*(u_t+k_ta-z_star_t-z_psi_t));
Eq_S13 = -y_t + c_y*c_t + i_y*i_t + z_g_t;

Eq_S14 = -k_t + (1-delta-r_n*psi/pii)/z_star/psi*u_t+ ...
         (1-delta)/z_star/psi*(k_ta-z_star_t-z_psi_t)+ ...
         (1-(1-delta)/z_star/psi)*(i_t+z_nu_t);

Eq_S15 = -q_t + inv_zeta*(i_t-i_ta+z_star_t+z_psi_t)-...
         inv_zeta*z_star*pii/r_n*(Et_i_t1-i_t+Et_z_star_t1+...
         Et_z_psi_t1)-z_nu_t+z_i_t;

Eq_S16 = -r_n_t + phi_r*r_n_ta+(1-phi_r)*...
        (0.25*phi_pi*(pi_t+pi_t_lg1+pi_t_lg2+pi_t_lg2a)+...
        phi_y*(y_t-y_ta))+z_r_t;

Eq_S17 = -Et_chi_t1 + r_E_t - Et_pi_t1;

%% Forecast Errors(9)
Eq_F1 = -c_t + Et_c_t1a + eta_c;
Eq_F2 = -lambda_t + Et_lambda_t1a + eta_lambda;
Eq_F3 = -w_t + Et_w_t1a + eta_w;
Eq_F4 = -chi_t + Et_chi_t1a + eta_chi;
Eq_F5 = -r_k_t + Et_r_k_t1a + eta_r_k;
Eq_F6 = -pi_t + Et_pi_t1a + eta_pi;
Eq_F7 = -i_t + Et_i_t1a + eta_i;
Eq_F8 = -z_b_t + Et_z_b_t1a + eta_z_b;
Eq_F9 = -z_star_t + Et_z_star_t1a + eta_z_star;

%% Persistent Shocks(11)
Eq_P1 = -z_b_t + rho_b*z_b_ta + epsilon_b;
Eq_P2 = -z_g_t + rho_g*z_g_ta + epsilon_g;
Eq_P3 = -z_w_t + rho_w*z_w_ta + epsilon_w;
Eq_P4 = -z_p_t + rho_p*z_p_ta + epsilon_p;
Eq_P5 = -z_nu_t + rho_nu*z_nu_ta + epsilon_nu;
Eq_P6 = -z_r_t + rho_r*z_r_ta + epsilon_r;
Eq_P7 = -z_z_t + rho_z*z_z_ta + epsilon_z;
Eq_P8 = -z_psi_t + rho_psi*z_psi_ta + epsilon_psi;
Eq_P9 = -z_i_t + rho_i*z_i_ta + epsilon_i;
Eq_P10 = -z_efp_t + rho_efp*z_efp_ta + epsilon_efp;
Eq_P11 = -z_nw_t + rho_nw*z_nw_ta + epsilon_nw;

%% Identities(3)
Eq_I1 = -z_star_t + z_z_t + alpha/(1-alpha)*z_psi_t;
Eq_I2 = -pi_t_lg1 + pi_ta;
Eq_I3 = -pi_t_lg2 + pi_t_lg1a;

%% n.eq = 40
System_of_Eq = [Eq_S1; Eq_S2; Eq_S3; Eq_S4; Eq_S5; ...
                Eq_S6; Eq_S7; Eq_S8; Eq_S9; Eq_S10; ...
                Eq_S11; Eq_S12; Eq_S13; Eq_S14; Eq_S15; ...
                Eq_S16; Eq_S17; ...
                Eq_F1; Eq_F2; Eq_F3; Eq_F4; Eq_F5; ...
                Eq_F6; Eq_F7; Eq_F8; Eq_F9;  ...
                Eq_P1; Eq_P2; Eq_P3; Eq_P4; Eq_P5; ...
                Eq_P6; Eq_P7; Eq_P8; Eq_P9; Eq_P10;  Eq_P11; ...
                Eq_I1; Eq_I2; Eq_I3];

npara = length(parameters);
neq  = length(System_of_Eq);  %  Num of stable and unstable Variables
nshock = length(Epsilon_t);
nend = length(Eta_t);         %  Num of Unstable Variables

C = zeros(neq,1);
GAM0j = -jacobian(System_of_Eq, Y_t);
GAM1j = jacobian(System_of_Eq, Y_ta);
PSI0j = jacobian(System_of_Eq, Epsilon_t);
PPIj  = jacobian(System_of_Eq, Eta_t);

%--------------------------------------
% output system matrix in fortran form
%--------------------------------------

% fortran(GAM0j);
% fortran(GAM1j);
% fortran(PSI0j);
% fortran(PPIj);

% z_star_t = z_z_t + alpha/(1-alpha)/z_psi_t
% z_g_t_tilde = z_g_t/(1-c_y-i_y)
% z_nw_t_tilde = z_nw_t/(1-z_star/r_E)

% c_y:   sample mean
% i_y:   sample mean
% l_bar:   sample mean
% r_n_bar: sample mean
% pi_bar: 1/4 for Japan

% z_star_bar : parameter
% psi_bar : parameter
% r_E_bar : parameter

% z_star = exp(0.01*z_star_bar);
% psi =exp(0.01*psi_bar);
% l = exp(0.01*l_bar);
% pii = exp(0.01*pi_bar);
% r_n = exp(0.01*r_n_bar);
% r_E = exp(0.01*r_E_bar);

% delta = 0.06
% alpha = 0.37
% lambda_w = 0.2
% lambda_i = 0.2

%---------------------------
% Solve the model by gensys
%---------------------------

addpath('./gensys');

% pos. mean of tab.2, p.28
sigma = 0.9114;
theta = 0.4907;
kai = 1.8996;
inv_zeta = 0.4113;
mu = 0.955;
phi_o_y = 0.9585;
gamma_w = 0.4309;
ksi_w = 0.5149;
gamma_p = 0.691;
ksi_p = 0.375;
phi_r = 0.453;
phi_pi = 1.754;
phi_y = 0.103;
z_star_bar = 0.289;
psi_bar = 0.407;
eta = 0.984;
n_k = 0.583;
mu_E = 0.0331*(1);
r_E_bar = 1.1845;

rho_b = 0.3442;
rho_g = 0.9245;
rho_w = 0.1625;
rho_p = 0.377;
rho_r = 0.2005;
rho_nu = 0.822;
rho_z = 0.2231;
rho_psi = 0.096;
rho_efp = 0.907;
rho_nw = 0.705;
rho_i = 0.976;

sigma_b = 1.46;
sigma_g = 1.01;
sigma_w = 0.61;
sigma_p = 0.58;
sigma_r = 0.61;
sigma_nu = 0.52;
sigma_z = 0.63;
sigma_psi = 0.59;
sigma_efp = 0.4135;
sigma_nw = 1.011;

delta = 0.06;
alpha = 0.37;
lambda_w = 0.2;
lambda_i = 0.2;

c_y = 1 - alpha;
i_y = alpha;
% l_bar = 1.0;
r_n_bar= 1.0;
pi_bar = 0.25;

z_star = exp(0.01*z_star_bar);
psi =exp(0.01*psi_bar);
% l = exp(0.01*l_bar);
pii = exp(0.01*pi_bar);
r_n = exp(0.01*r_n_bar);
r_E = exp(0.01*r_E_bar);

GAM0 = eval(GAM0j);
GAM1 = eval(GAM1j);
PSI0 = eval(PSI0j);
PPI  = eval(PPIj);

[T1,TC,T0,fmat,fwt,ywt,gev,RC,loose] = gensys(GAM0,GAM1,C,PSI0,PPI,1);

%-------------------
% Impulse Response
%-------------------
%2014.8.15�����C��Z0
Z0 = [z_star_t + y_t - y_ta;
      z_star_t + c_t - c_ta;
      z_star_t + z_psi_t + i_t - i_ta;
      z_star_t + w_t - w_ta;
      l_t;
      pi_t;
      -z_psi_t + z_nu_t - z_nu_ta;
      r_n_t;
      Et_chi_t1;
      z_star_t + b_t - b_ta;
      k_t - k_ta + z_star_t + z_psi_t;
      n_t - n_ta + z_star_t;
      q_t - q_ta - z_psi_t;
      Et_r_E_t1];

nvar = size(Z0, 1);
%nirf = 20; %2014.8.15�����C��
Z = jacobian(Z0, Y_t); 
Zlg = jacobian(Z0, Y_ta); 
ZZ = [Z,Zlg];

ZZ1 = eval(ZZ);

sig_chol = eye(nshock); 

% titlestr = {'Zb :', 'Zg :', 'Zw :', 'Zp :', 'Z_\nu :',...
%             'Zr :', 'Zz :', 'Z_\phi :', 'Z_i:', 'Zefp :', 'Znw :'};
% 
% ystr = {'Yt', 'Ct', 'It', 'Wt', 'Lt', 'Pt', 'Pi/Pt',...
%         'rn_t', 'Bt','Kt','Nt','Qt', 'R^E_t'};
    
    
  yyirf_total  = zeros(nirf,nvar,nshock);  % 2014/08/10 �ǉ�

for sh_ind = 1:nshock
    
%  figure('Name', 'Impulse Responses of Structual Shock' ); % 2014/08/10 �C��
% figure('Name',strcat(titlestr(sh_ind)) );
    
    impact = sig_chol(:,sh_ind);
    yyirf  = zeros(nirf,nvar);
    dyyirf  = zeros(nirf,nvar);
    ss = T0*impact;
    s  = [ss;zeros(size(ss,1),1)];
    dyyirf(1,:) = ZZ*s;
    yyirf(1,:) = dyyirf(1,:);
    
    for t = 2:nirf
        ss1 = T1*ss;     
        s = [ss1;ss];
        dyyirf(t,:) = (ZZ*s)';
        yyirf(t,1:4) = yyirf(t-1,1:4)+dyyirf(t,1:4);
        yyirf(t,5:nvar) = dyyirf(t,5:nvar);%2014.8.15�����C��
        ss = ss1;
    end
%     for j = 1:9
%       subplot(3, 3, j)
%       plot(1:nirf, yyirf(:,j),'b')
%       title(strcat(titlestr(sh_ind),ystr(j)))
%     end
    
    yyirf_total(:,:,sh_ind)  = yyirf;  % 2014/08/10 �ǉ�
    
%   display('press any key to change to next graph');  
%   w = waitforbuttonpress;
end

