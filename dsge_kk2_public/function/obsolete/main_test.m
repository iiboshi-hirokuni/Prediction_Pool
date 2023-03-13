%=============================
% Main Program
%=============================

cd 'H:/work/Hprograming/matlab/DSGE_SWM'
addpath('./function');
addpath('./gensys');
addpath('./csmin');

%---------------------
% Path Specification
%---------------------
priopath = './prior/';
postpath = './posterior/';
datapath = './data/';
resupath = './results/';

%---------------------
% Model Specification
%---------------------

nvar   = 7;  % num of  observed variables
nshock = 9;
npara  = 27;
ngrid  = 100; % number of grid points for density plots

runname = 'm1011'
datafilename = 'm1sim1_test.csv';
runpath = runname;

trspecfile = strcat(priopath,'_trspec.csv');
trspec = csvread(trspecfile,1,1);
%    Each row has the following specification:
% 
%     tr~a~b~c
% 
%     tr parameter transformation type
%         0: no transformation needed
%         1: [a,b] -> [-1,1] -> [-inf,inf] by (1/c)*c*z/sqrt(1-c*z^2)
%         2: [0,inf] -> [-inf,inf] by b + (1/c)*ln(para[i]-a);
%     a  transformation argument a (usually lower bound)
%     b  transformation argument b (usually upper bound)
%     c  transformation argument c

para_names={'h','sigma_c','sigma_L','phi','phi0','psi','gam_p',...
            'gam_w','xi_p','xi_w','rho_m','mu_pi','mu_y','rho_z',...
            'rho_c','rho_g','rho_L','rho_i','e_c','e_inv','e_q',...
            'e_L','e_w','e_z','e_p','e_g','e_m'};

%-------------------
% loading DATA
%-------------------

series_YT = csvread(strcat(datapath, datafilename), 1, 1);

nobs = size(series_YT, 1);  % number of observations 
yy_m = mean(series_YT, 1);
for i = 1:nobs
  YY(i,:) = series_YT(i,:) - yy_m;
end

%dataplot

%-------------------------------------
% Import the information for the prior density
%-------------------------------------

priorfile = strcat(priopath, runname,'prt.csv');
prior = csvread(priorfile, 1, 1);

pmean  = prior(:,1);
pstdd  = prior(:,2);
pshape = prior(:,3);
pmask  = prior(:,4);
pfix   = prior(:,5);
pmaskinv = 1-pmask;
pshape   = pshape.*pmaskinv;

varargin = struct('pmaskinv',pmaskinv,'pfix',pfix,'pmask',pmask,'pmean',...
     pmean,'pstdd',pstdd,'pshape',pshape,'data',YY,'trspec',trspec,'nshock',nshock);

% block size for MH sampling
nblock = 11;
% number of simulation per block
nsim   = 100;

nburnin_block = 1;
nthin = 10;

%-------------------
% test
%--------------------

h       = 0.5;
sigma_c = 1;
sigma_L = 1;
phi     = 4;
phi0    = 2;
psi     = 0.15;
gam_p   = 0.3;
gam_w   = 0.7;
xi_p    = 0.88; 
xi_w    = 0.62;
rho_m   = 0.8;
mu_pi   = 1.7;
mu_y    = 0.12;
rho_z  =  0.6;
rho_c   =  0.5;
rho_g   =  0.9;
rho_L   =  0.9; 
rho_i   =  0.75;
e_c     =  2;
e_inv   = 1.8;
e_q     = 6;
e_L     = 1;
e_w     = 10;
e_z     = 1;
e_p     = 10;
e_g     = 1;
e_m     = 0.2;

para    = [ h ; sigma_c ; sigma_L ; phi ; phi0 ; psi ; gam_p ;
            gam_w ; xi_p ;  xi_w ; rho_m ;
            mu_pi ; mu_y ; rho_z ; rho_c ; rho_g ; rho_L ; rho_i ;
            e_c ; e_inv ; e_q ;  e_L ;  e_w ; e_z ; e_p ; e_g ; e_m ];

para_sv = para

%******************************************************** 
% Calculate posterior density at starting values

% test of Kalman filter

fprintf('Prior Density at Starting Value:\n');
priodens(para,pmean,pstdd,pshape)

[lnpY,retcode,obsmean,obsvar,shock] = evaldsge(para,YY,nshock);

fprintf('Posterior at Starting Value: %f\n', lnpY);


% test of CEEpm.m

nac = nvar + nvar*(nvar+1)/2 + nvar^2;
nvd = nvar * 9;

nn = 4;   % numbers of Horizon

        
        [GAM0_c, GAM0_i, GAM0_q, GAM0_l, GAM0_w, GAM0_z, GAM0_p, GAM0_g, GAM0_m, ...
          GAM1_c, GAM1_i, GAM1_q, GAM1_l, GAM1_w, GAM1_z, GAM1_p, GAM1_g, GAM1_m ...
            ] = nkbcmom(para, nn);
        
        GAM0 = GAM0_c + GAM0_i + GAM0_q + GAM0_l ...
            + GAM0_w + GAM0_z + GAM0_p + GAM0_g + GAM0_m;
        GAM1 = GAM1_c + GAM1_i + GAM1_q + GAM1_l ...
            + GAM1_w + GAM1_z + GAM1_p + GAM1_g + GAM1_m;
        
        % Variance Decompositions
        vd  = [(diag(GAM0_c)./diag(GAM0)), (diag(GAM0_i)./diag(GAM0)),...
               (diag(GAM0_q)./diag(GAM0)), (diag(GAM0_l)./diag(GAM0)),...
               (diag(GAM0_w)./diag(GAM0)), (diag(GAM0_z)./diag(GAM0)),...
               (diag(GAM0_p)./diag(GAM0)), (diag(GAM0_g)./diag(GAM0)),...
               (diag(GAM0_m)./diag(GAM0))];
        %vdsim[j,.] =  vec(vd')';
        %hvddraw(j,:) = reshape(vd', nvd, 1);
hvddraw = reshape(vd', nvd, 1)'
        
        % Autocovariances
        corr0  = corrvc(GAM0);
        stddev = sqrt(diag(GAM0));
        corr1  = GAM1 ./ (stddev*stddev');
        
        %hacdraw(j,:) = [stddev' , vech(corr0)' ,reshape(corr1',nvar*nvar, 1)'];
hacdraw = [stddev' , vech(corr0)' ,reshape(corr1',nvar*nvar, 1)']


%-------------------------------------
% Estimation of DSGE model using data
%-------------------------------------

%% Maximize Posterior Density
%  Maximize the posterior density, calculate Hessian, display results

%CEEpm

%% Compute Hessian at Posterior Mode
%  compute square root inverse hessian

%CEEhess
%sqrt(diag(inv(hessian)))

%% Metropolis-Hastings RW Simulation
%  Metropolis-Hastings step

% CEEmh

% Variance Decomposition

% CEE_ac_n

% Inpulse Response Functions

% CEE_irf

