%=====================================================================
%  Main Program
%  DSGE model with Finicial Friction
%
%=====================================================================

clear all
clc

warning('off','all'); % 警告の非表示

global ZZ YY ZZ_leverage;

addpath('./function');
addpath('./function_KK');
addpath('./gensys');
addpath('./csmin');

%---------------------------------------------------------------------
% Path Specification
%---------------------------------------------------------------------
priopath = './prior/';
postpath = './posterior/';
datapath = './data/';
resupath = './results/';

%---------------------------------------------------------------------
% Model Specification
%---------------------------------------------------------------------

nvar   = 10;  % num of  observed variables
nshock = 11;
npara  = 41;
ngrid  = 100; % number of grid points for density plots

%runname = 'm1011'
%datafilename = 'm1sim1.csv';
runname = 'FF2011'
datafilename = 'kaihatsu_kurozumi_data_1.csv';

runpath = runname;
trspecfile = strcat(priopath,'_trspec.csv');
trspec = csvread(trspecfile,1,1);

para_names = {'sigma','theta','kai','inv_zeta',...
     'mu','phi_o_y','gamma_w','ksi_w', ...
     'gamma_p','ksi_p','phi_r','phi_pi','phi_y','z_star_bar','psi_bar', ...
     'eta','n_k','mu_E','r_E_bar', ...
     'rho_b','rho_g','rho_w','rho_p','rho_r','rho_nu','rho_z','rho_psi', 'rho_i',...
     'rho_efp','rho_nw', ...
     'sigma_b','sigma_g','sigma_w','sigma_p','sigma_r', ...
     'sigma_nu','sigma_z','sigma_psi','sigma_i','sigma_efp','sigma_nw'};

%---------------------------------------------------------------------
% loading DATA
%---------------------------------------------------------------------

series_YT = csvread(strcat(datapath, datafilename), 1, 1);

nobs = size(series_YT, 1);  % number of observations 
yy_m = mean(series_YT, 1);


YY = series_YT;

%---------------
% plotting data
%---------------

%  dataplot
%w = waitforbuttonpress;

%---------------------------------------------------------------------
% Import the information for the prior density
%---------------------------------------------------------------------

priorfile = strcat(priopath, runname,'prt.csv');
prior = csvread(priorfile, 1, 1);

pmean  = prior(:,1);
pstdd  = prior(:,2);
pshape = prior(:,3);
pmask  = prior(:,4);
pfix   = prior(:,5);
pmaskinv = 1-pmask;
pshape   = pshape.*pmaskinv;

%% 観測方程式の行列の設定

[ ZZ, ZZ_leverage] = make_zz();

varargin = struct('pmaskinv',pmaskinv,'pfix',pfix,'pmask',pmask,...
     'pmean',pmean,'pstdd',pstdd,'pshape',pshape,'data',YY,...
     'trspec',trspec,'nshock',nshock, 'ZZ', ZZ);

% block size for MH sampling
nblock = 50;
% number of simulation per block
nsim   = 1000;

nburnin_block = round(0.5*nblock);   % blocks for burn-in
nthin = 1;                           % number of thining
 cc = 0.08;                          % adjustment coefficient of random walk MH algorithm
                                     % to set cc to become acceptance rate 25% (or rejection rate 75%)

%=====================================================================
% Estimation of DSGE model using data
%=====================================================================

%---------------------------------------------------------------------
% Generate valid draws from the prior of the DSGE model
%---------------------------------------------------------------------
%
%   DSGE_pri
% 
%  plot_type = 'prior'
%  plot_dist % 事前後分布の作成
%
%---------------------------------------------------------------------
% Maximize Posterior Density
% Maximize the posterior density, calculate Hessian, display results
%---------------------------------------------------------------------

%   DSGE_pm

%---------------------------------------------------------------------
% Compute Hessian at Posterior Mode
% compute square root inverse hessian
%---------------------------------------------------------------------

%  DSGE_hess

%---------------------------------------------------------------------
% Metropolis-Hastings RW Simulation
% Metropolis-Hastings step
%---------------------------------------------------------------------

 use_post_mode = 1;  % Yes-> 1, No->0
 use_post_hess = 0;
 
%               DSGE_mh
% % 
          stats_sample_para
% %     
    plot_type = 'posterior'
%     plot_dist  %　事後分布の作成
%     plot_trace   %　事後分布のサンプリングのトレースの作成

%    Marginal Likelihood
%     mcmcdd

%---------------------------------------------------------------------
% Impulse Response Functions
%---------------------------------------------------------------------

%       DSGE_irf
% 
% %---------------------------------------------------------------------
% % Forecasting observable variable
% %---------------------------------------------------------------------
%      
%              DSGE_Forecast
        
% %---------------------------------------------------------------------
% %  Calculating Prediction Density
% %---------------------------------------------------------------------
%      
%              KK_Pred_Pool
               

% %---------------------------------------------------------------------
% %  Historical Decomposition
% %---------------------------------------------------------------------
    
% for i = 1:nobs
%  YY(i,:) = series_YT(i,:) - yy_m;
% end
% 
%     DSGE_hist_decomp


     
     