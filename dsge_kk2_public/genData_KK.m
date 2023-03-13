
%*************************************************************************
%  Simulate the DSGE model 
%*************************************************************************

%---------------------
% Path Specification
%---------------------

addpath('./function');
addpath('./function2');
addpath('./gensys');

priopath = './prior/';
datapath = './data/';

outputfilename = 'm1sim1.csv';

%---------------------
% Model Specification
%---------------------

nvar   = 10;
nstate = 40;

Tsimtotal = 200;
Tsimburn  = 200;
%Tsim = [200 ; 500];

% Load model parameters
lppara = csvread(strcat(priopath,'msim_kk_par.csv'), 1, 1);

para = lppara(:,1); 

sigma      = para(1,1);
theta      = para(2,1);
kai        = para(3,1);
inv_zeta   = para(4,1);
mu         = para(5,1);
phi_o_y    = para(6,1);
gamma_w    = para(7,1);
ksi_w      = para(8,1);
gamma_p    = para(9,1);
ksi_p      = para(10,1);
phi_r      = para(11,1);
phi_pi     = para(12,1);
phi_y      = para(13,1);
z_star_bar = para(14,1);
psi_bar    = para(15,1);
eta        = para(16,1);
n_k        = para(17,1);
mu_E       = para(18,1);
r_E_bar    = para(19,1);

sigma_b    = para(30,1);
sigma_g    = para(31,1);
sigma_w    = para(32,1);
sigma_p    = para(33,1);
sigma_r    = para(34,1);
sigma_nu   = para(35,1);
sigma_z    = para(36,1);
sigma_psi  = para(37,1);
sigma_efp  = para(38,1);
sigma_nw   = para(39,1);

l_bar   = 1.0;
r_n_bar = 1.0;
pi_bar  = 0.25;

% Solve the DSGE model */
%
% retcode = -1 : non existence                
%         = 0  : existence and uniqueness     
%         = 1  : existence and non-uniqueness

[T1,TC,T0,RC] = dsgesolv(para);      %T1; TC; T0; RC;

% Simulate DSGE model
% Composition of xt: [x(t),pi(t),Et(x(t+1)),Et(pi(t+1)),m(t),deltaM(t),z(t)]'

%ZZ = zeros(9, size(T1, 1)*2);
%ZZ(1,14) = 1; ZZ(1,29)=1; ZZ(1,54) = -1;                % output gap  1
%ZZ(2,1)  = 1; ZZ(2,29)=1; ZZ(2,41) = -1;                % consumption 2
%ZZ(3,15) = 1; ZZ(3,26)=1; ZZ(3,29) =  1; ZZ(3,55) = -1; % investment  3
%ZZ(4,3)  = 1; ZZ(4,29)=1; ZZ(4,43) = -1;                % wage        4
%ZZ(5,30) = 1;                                           % labor       5
%ZZ(6,11) = 4;                                           % inflation   6
%ZZ(7,23) = 4; ZZ(7,26)=-4; ZZ(7,63) = -4;               % investment deflator 7
%ZZ(8,17) = 4;                                           % nominal rate   8
%ZZ(9,4)  = 1; ZZ(9,29)=1;  ZZ(9,44) = -1;               % real borrowing 9

ZZ = make_zz();
DD = zeros(nvar, 1);
DD(1, 1) = z_star_bar;
DD(2, 1) = z_star_bar;
DD(3, 1) = z_star_bar+psi_bar;
DD(4, 1) = z_star_bar;
DD(5, 1) = l_bar;
DD(6, 1) = pi_bar;
DD(7, 1) = -psi_bar;
DD(8, 1) = r_n_bar;
DD(9, 1) = 0;
DD(10, 1) = z_star_bar;

% Initializations 
%
yT = zeros(Tsimtotal, nvar);
xt = zeros(nstate,1);

for t = 1:Tsimburn+Tsimtotal
   % calculate realization of epst
   epst  = [randn()*sigma_b;
            randn()*sigma_g;
            randn()*sigma_w;
            randn()*sigma_p;
            randn()*sigma_nu;
            randn()*sigma_r;
            randn()*sigma_z;
            randn()*sigma_psi;
            randn()*sigma_efp;
            randn()*sigma_nw];
   xtlg = xt;
   xt = T1*xt + TC + T0*epst;
   % Update vector x(t)
   if t > Tsimburn
      yT(t-Tsimburn,:) = DD' + (ZZ*[xt;xtlg])';
   end
end

%yT = yT(Tsimburn+1:Tsimtotal+Tsimburn, :);

fprintf('Data Simulated Successfully\n')
fprintf('Parameters (theta):\n')
para'
fprintf('Mean of YT:\n');
mean(yT, 1)
fprintf('StD of YT:\n');
std(yT, 0, 1)
fprintf('RC=\n');
RC
fprintf('1: success  0: fail\n');

%---------------------
% Save simulated data
%---------------------

osimdata = strcat(datapath, outputfilename);

rowlabel = [1:size(yT,1)]';
collabel = {'Time', 'Output', 'Consumption','Investment', 'Wage', 'Labor',  'Inflation',...
            'Investment_deflator', 'Nominal_rate', 'Output_Gap', 'Real_borrowing' };

util_csvwrite(osimdata, [rowlabel, yT], collabel);
fprintf('Data saved\n');

