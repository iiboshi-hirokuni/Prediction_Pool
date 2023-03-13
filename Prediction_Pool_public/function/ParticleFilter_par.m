function [ re_lambda, log_lik,mean_lambda ]=...
    ParticleFilter_par(PredDen_1, PredDen_2, T0_Forecast,...
    h_Forecast, particles, rho, mu, sigma, shock_save,ncores) 

 pd              =   'Normal';

% particles = nsim;
npar = particles/ncores;
t       = size(PredDen_1,1);                 % Number of period
xt      = zeros(t,2,particles);                  % Particles of factors (3*1) at period t
% lambda  = zeros(t, particles);
% et      = zeros(t, particles);                  % Shock e_t ~ N(0,1) coresponding to each particles h_t
% wt0      = zeros(t, particles);                  % Weight of each particles
wt      = zeros(npar,ncores);
xt_temp = zeros(2,particles);
      
  PredDen = [ PredDen_1  PredDen_2];
  P.rho = rho; 
  P.mu  = mu;
  P.sigma = sigma;

 shock = squeeze(shock_save(1,:,:));
 
 xt_Resamp = zeros(2,particles);
 xt_Resamp(1,:) = (1-rho)*mu+rho*rand(particles,1) + sqrt(1 - rho^2)*sigma*shock(:,1);
  

 parfor k =1:ncores                          
          [wt(:,k), a]= sub_parallel_2model(P, PredDen, shock(1+(k-1)*npar:k*npar),...
                                              xt_Resamp(:,1+(k-1)*npar:k*npar),npar,T0_Forecast,2);         
          x_temp(:,:,k)     =  a;      %  sampling of x_t
 end

  for k = 1:ncores
         prob_x(1,1+(k-1)*npar:k*npar)  = wt(:,k);  
         xt_temp(:,1+(k-1)*npar:k*npar) = squeeze(x_temp(:,:,k)); 
  end

    lik_t = mean(prob_x);          
    log_lik = log(lik_t); 
  
%    parfor k = 1:ncores
%       xt_Resamp1(:,:,k) = resample1(xt_temp, prob_x, particles,ncores);
%    end
%    
%    for k = 1:ncores
%       xt_Resamp(:,1+(k-1)*npar:k*npar) =  squeeze(xt_Resamp1(:,:,k)) ;
%    end    
   
%    
%    for j=1:particles
%         xt_Resamp(:,j) = mnrnd(1, prob_x/sum(prob_x)) * xt_temp';
%     end
    xt_Resamp = resample(xt_temp, prob_x/sum(prob_x), particles);
    xt(T0_Forecast,:,:) = xt_Resamp;

for t = T0_Forecast+1:1:t-h_Forecast           
      shock = squeeze(shock_save(t,:,:));
      
    parfor k =1:ncores                          
        [wt(:,k), a]= sub_parallel_2model(P, PredDen, shock(1+(k-1)*npar:k*npar),...
                                          xt_Resamp(:,1+(k-1)*npar:k*npar),npar,t,2);         
         x_temp(:,:,k)     =  a;      %  sampling of x_t
    end    
  
   for k = 1:ncores
         prob_x(1,1+(k-1)*npar:k*npar)  = wt(:,k);  
         xt_temp(:,1+(k-1)*npar:k*npar)    = squeeze(x_temp(:,:,k)); 
    end  
    

%     for j=1:particles             
%         xt_Resamp(:,j) = mnrnd(1, prob_x/sum(prob_x)) * xt_temp';
%     end
    
%     parfor k = 1:ncores
%       xt_Resamp1(:,:,k) = resample1(xt_temp, prob_x, particles,ncores);
%    end
%    
%    for k = 1:ncores
%       xt_Resamp(:,1+(k-1)*npar:k*npar) =  squeeze(xt_Resamp1(:,:,k)) ;
%    end    

    xt_Resamp = resample(xt_temp, prob_x/sum(prob_x), particles);
    xt(t,:,:) = xt_Resamp;
    
    % Step 4: Calculating likelihood
    lik_t = mean(prob_x);          
    log_lik = log_lik + log(lik_t);  % likelihood of full period
end          
  
% Sampling one partilce from the distribution   
re_lambda = squeeze(xt(:,2,:));

mean_lambda = mean(squeeze(xt(:,2,:)),2);
  
 
  
  

  
  
  
  
  